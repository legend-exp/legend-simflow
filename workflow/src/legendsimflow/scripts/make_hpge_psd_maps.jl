# Copyright (C) 2025 Giovanna Saleh <giovanna.saleh@phd.unipd.it>,
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

# grid spacing in mm
const DEFAULT_GRID_SIZE_IN_MM = 5.0
# crystal axis angles in degrees (<001> and <110>)
const CRYSTAL_AXIS_ANGLES = [0, 45]
# SSD adaptive-mesh refinement thresholds as fractions of the crystal radius
# matches current SSD behaviour
const DEFAULT_REFINEMENT_LIMITS = [0.2, 0.1, 0.05, 0.02]
# nr of pixels for padding around the map to avoid grid edge effects
const DEFAULT_PADDING = 3

using LegendHDF5IO
using ArgParse
using PropDicts
using Printf
using Unitful
using LegendSimflow

"""
    main()

Generate HPGe ideal pulse-shape libraries or drift-time maps for a detector
and save them to an LH5 file.
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
        "--info-file"
        help = "Path to output YAML file with the SSD-modeling provenance scalars (optional)"
        default = nothing
    end
    @add_arg_table s begin
        "--opv"
        help = "detector operational voltage in V (defaults to metadata value)"
    end
    @add_arg_table s begin
        "--grid-size-in-mm"
        help = "Spacing of the (r, z) grid, in mm"
        arg_type = Float64
        default = DEFAULT_GRID_SIZE_IN_MM
    end
    @add_arg_table s begin
        "--refinement-limits"
        help = "SSD adaptive-mesh refinement thresholds, one per refinement pass"
        arg_type = Float64
        nargs = '+'
        default = DEFAULT_REFINEMENT_LIMITS
    end
    @add_arg_table s begin
        "--padding"
        help = "Number of pixel layers added around the (r, z) grid"
        arg_type = Int
        default = DEFAULT_PADDING
    end
    @add_arg_table s begin
        "--compute-drift-time"
        help = "Flag to compute drift time maps instead of ideal pulse shapes (optional; default: false)"
        action = :store_true
    end
    @add_arg_table s begin
        "--impurity"
        help = "YAML file with the tuned impurity-profile slope and depletion voltage, keyed by detector (optional; metadata values used if absent)"
        default = nothing
    end
    parsed_args = parse_args(s)

    det = parsed_args["detector"]
    meta_path = parsed_args["metadata"]
    output_file = parsed_args["output-file"]
    info_file = parsed_args["info-file"]

    !isnothing(info_file) && isfile(info_file) && error("Info file already exists")
    isfile(output_file) && error("Output file already exists")


    raw_opv = parsed_args["opv"]
    opv_val = isnothing(raw_opv) ? nothing : parse(Float32, raw_opv)
    meta, xtal, opv_val = load_detector_metadata(meta_path, det, opv_val)


    # tuned impurity profile: rescale its non-constant part by `slope` and
    # match the simulated depletion voltage (in V) to `depletion_voltage`
    vdep = nothing
    impurity_file = parsed_args["impurity"]
    if !isnothing(impurity_file)
        impurity = readprops(impurity_file)[Symbol(det)]
        vdep = impurity.depletion_voltage
        @info "Adjusting impurity profile parameters with slope = $(impurity.slope)"
        xtal.impurity_curve.parameters = adjust_impurity_pars(xtal.impurity_curve.parameters, impurity.slope)
    end

    grid_size = parsed_args["grid-size-in-mm"] / 1000.0
    ref_limits = parsed_args["refinement-limits"]
    padding = parsed_args["padding"]

    @info "using ref limits $ref_limits"
    sim, info = setup_hpge_simulation(meta_path, meta, xtal, opv_val, T, ref_limits, vdep = vdep)

    output = nothing
    for a in CRYSTAL_AXIS_ANGLES

        if !parsed_args["compute-drift-time"]
            result = compute_ideal_pulse_shape_lib(sim, meta, T, a, false, grid_size, padding)

            key = Symbol("waveform_$(lpad(string(a), 3, '0'))_deg")
        else
            result = compute_drift_time_map(sim, meta, T, a, grid_size, padding)

            key = Symbol("drift_time_$(lpad(string(a), 3, '0'))_deg")
        end

        if output === nothing
            output = Dict{Symbol,Any}(pairs(result))
        else
            output[key] = result[key]
        end
    end


    @info "Saving to disk..."
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        @info "Creating output directory: $output_dir"
        mkpath(output_dir)
    end

    lh5open(output_file, "cw") do f
        return f[det] = (; output...)
    end

    if !isnothing(info_file)
        @info "Saving SSD-modeling info to $info_file..."
        writeprops(info_file, PropDict(info))
    end
end


main()
