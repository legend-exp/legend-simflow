# Copyright (C) 2026 Giovanna Saleh <giovanna.saleh@phd.unipd.it>,
#                    Luigi Pertoldi <gipert@pm.me>,
#                    Toby Dixon <toby.dixon.23@ucl.ac.uk> and
#                    David Hervas <david.hervas@tum.de>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# grid spacing in meters
const DEFAULT_GRID_SIZE = 0.005
# crystal axis angles in degrees (<001> and <110>)
const CRYSTAL_AXIS_ANGLES = [0, 45]
# SSD adaptive-mesh refinement thresholds as fractions of the crystal radius
# matches current SSD behaviour
const DEFAULT_REFINEMENT_LIMITS = [0.2, 0.1, 0.05, 0.02]
# nr of pixels for padding around the map to avoid grid edge effects (default; can be overridden via metadata settings file)
const DEFAULT_PADDING = 3

# scan grid, as ranges of depletion-voltage shift in V and of dimensionless
# impurity-profile scaling factor; overridable via the scan settings file
const DEFAULT_DEPV_SHIFTS = -1000.0:20.0:-20.0
const DEFAULT_SLOPES = -0.9:0.2:3.0

using LegendHDF5IO
using ArgParse
using PropDicts
using Printf
using Unitful
using LegendSimflow
using SolidStateDetectors

"""Parse a `"start:step:stop"` scan range from the settings file."""
function parse_range(str::AbstractString)
    parts = split(str, ':')
    length(parts) == 3 || error("expected a \"start:step:stop\" range, got \"$str\"")
    vals = parse.(Float64, parts)
    return vals[1]:vals[2]:vals[3]
end

"""
    main()

Scan the ideal HPGe pulse shape library over the impurity profile and the
depletion voltage, and write it to an LH5 file.

For each scaling factor of the impurity profile (the `slope` of the scan) the
crystal impurity parameters are rescaled with [`adjust_impurity_pars`](@ref)
and the detector is simulated again. For each depletion voltage of the scan the
impurities are then rescaled to match it, and the waveform map is computed at
every crystal axis angle.

# Inputs

Command line arguments:

- `--detector`: HPGe detector name, e.g. `V05261B`
- `--metadata`: path to legend-metadata, source of the detector and crystal
  parameters
- `--opv`: operational voltage in V; read from the metadata when absent
- `--ssd-settings`: YAML file with the SSD simulation settings
  (`grid_size_in_mm`, `ssd_refinement_limits`, `padding`), see
  [`setup_hpge_simulation`](@ref). Built-in defaults are used when the file is
  absent
- `--scan-settings`: YAML file with the scan grid, given as two Julia ranges
  written `"start:step:stop"`: `depv_shift`, the depletion voltage relative to
  the operational voltage in V, and `slope`, dimensionless. Built-in defaults
  are used when the file is absent
- `--output-file`: path of the output LH5 file, which must not exist yet

# Output

A single LH5 file holding one group named after the detector:

```
<detector>/
├── psl_scan/
│   ├── slope_1/                  # one group per impurity slope, in scan order
│   │   ├── dep_1/                # one group per depletion voltage
│   │   │   ├── r                 # radial axis, in m
│   │   │   ├── z                 # axial axis, in m
│   │   │   ├── dt                # waveform sampling period, 1 ns
│   │   │   ├── waveform_000_deg  # normalized to 1, see below for the shape
│   │   │   └── waveform_045_deg
│   │   └── dep_2/ ...
│   └── slope_2/ ...
└── info/
    ├── slope_min                 # first slope of the scan
    ├── slope_step
    ├── dep_min                   # first depletion voltage of the scan, in V
    └── dep_step                  # in V
```

The groups are numbered from 1 in scan order, so the values behind `slope_i`
and `dep_j` are `slope_min + (i - 1) * slope_step` and
`dep_min + (j - 1) * dep_step`.

The waveform arrays are indexed `[time, z, r]` in Julia. Julia writes them in
column-major order, so a row-major reader such as `h5py` or `lgdo` sees the
reversed shape `(n_r, n_z, n_time)`.
"""
function main()
    T = Float64

    s = ArgParseSettings()

    @add_arg_table s begin
        "--detector"
        help = "HPGe detector name"
        required = true
    end
    @add_arg_table s begin
        "--metadata"
        help = "Path to legend-metadata"
        required = true
    end
    @add_arg_table s begin
        "--output-file"
        help = "Path to output LH5 file"
        required = true
    end
    @add_arg_table s begin
        "--opv"
        help = "detector operational voltage in V (defaults to metadata value)"
    end
    @add_arg_table s begin
        "--ssd-settings"
        help = "Path to ssd settings YAML file (optional; built-in defaults used if absent or missing)"
        default = nothing
    end
    @add_arg_table s begin
        "--scan-settings"
        help = "Path to scan settings YAML file (optional; built-in defaults used if absent or missing)"
        default = nothing
    end
    parsed_args = parse_args(s)

    det = parsed_args["detector"]
    meta_path = parsed_args["metadata"]
    output_file = parsed_args["output-file"]

    isfile(output_file) && error("Output file already exists")

    # extract the metadata
    raw_opv = parsed_args["opv"]
    opv_val = isnothing(raw_opv) ? nothing : parse(T, raw_opv)

    meta, xtal, opv_val = load_detector_metadata(meta_path, det, opv_val)

    # Load optional simulation settings, falling back to built-in defaults.
    # The settings file path is passed via --ssd-settings and applies globally
    # to all detectors and voltages.
    ssd_settings = parsed_args["ssd-settings"]
    sim_cfg = (!isnothing(ssd_settings) && isfile(ssd_settings)) ? readprops(ssd_settings) : PropDict()
    grid_size = get(sim_cfg, :grid_size_in_mm, DEFAULT_GRID_SIZE * 1000) / 1000.0
    ref_limits = get(sim_cfg, :ssd_refinement_limits, DEFAULT_REFINEMENT_LIMITS)
    padding = get(sim_cfg, :padding, DEFAULT_PADDING)

    # extract the settings of the scan
    scan_settings = parsed_args["scan-settings"]
    scan_cfg = (!isnothing(scan_settings) && isfile(scan_settings)) ? readprops(scan_settings) : PropDict()

    depv_shifts =
        haskey(scan_cfg, :depv_shift) ? parse_range(scan_cfg.depv_shift) : DEFAULT_DEPV_SHIFTS
    slopes = haskey(scan_cfg, :slope) ? parse_range(scan_cfg.slope) : DEFAULT_SLOPES


    # loop over slopes
    output = Dict{Symbol,Any}()
    base_xtal = deepcopy(xtal)

    time_setup = 0
    time_rescale = 0

    time_drift = 0

    for (sidx, slope) in enumerate(slopes)

        t0 = time()
        xtal.impurity_curve.parameters = adjust_impurity_pars(base_xtal.impurity_curve.parameters, slope)

        sim, _ =
            setup_hpge_simulation(meta_path, meta, xtal, opv_val, T, ref_limits, vdep = opv_val + first(depv_shifts))

        time_setup += time() - t0

        output[Symbol("slope_$sidx")] = Dict{Symbol,Any}()

        for (didx, depv_shift) in enumerate(depv_shifts)
            depv = opv_val + depv_shift

            t0 = time()
            scale = adjust_impurity_and_electric_potential_to_match_depletion!(sim,
                depv,
                check_for_depletion = false,
                reconverge_electric_potential = false
            )

            calculate_electric_field!(sim)
            time_rescale += time() - t0

            output[Symbol("slope_$sidx")][Symbol("dep_$didx")] = nothing

            t0 = time()
            for a in CRYSTAL_AXIS_ANGLES
                result = compute_ideal_pulse_shape_lib(sim, meta, T, a, false, grid_size, padding)

                key = Symbol("waveform_$(lpad(string(a), 3, '0'))_deg")
                if output[Symbol("slope_$sidx")][Symbol("dep_$(didx)")] === nothing
                    output[Symbol("slope_$sidx")][Symbol("dep_$(didx)")] = Dict{Symbol,Any}(pairs(result))
                else
                    output[Symbol("slope_$sidx")][Symbol("dep_$(didx)")][key] = result[key]
                end
            end

            time_drift += time() - t0
        end
    end

    @info "Timing summary:"
    @info "  Setup time: $(time_setup) seconds"
    @info "  Rescale time: $(time_rescale) seconds"
    @info "  Drift time: $(time_drift) seconds"

    t0 = time()
    @info "Saving to disk..."
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        @info "Creating output directory: $output_dir"
        mkpath(output_dir)
    end
    # reformat

    dict_to_namedtuple(d::AbstractDict) =
        (; (k => v isa AbstractDict ? dict_to_namedtuple(v) : v for (k, v) in d)...)


    output = dict_to_namedtuple(output)
    output = (
        psl_scan = output,
        info = (
            slope_min = first(slopes),
            slope_step = step(slopes),
            dep_min = opv_val + first(depv_shifts),
            dep_step = step(depv_shifts)
        )
    )

    lh5open(output_file, "cw") do f
        return f[det] = output
    end
    time_write = time() - t0
    @info "Saved to disk in $(time_write) seconds"

end


main()
