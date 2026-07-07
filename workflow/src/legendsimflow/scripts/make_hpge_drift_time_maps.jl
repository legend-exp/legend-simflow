# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

# grid spacing in meters (default; can be overridden via metadata settings file)
const DEFAULT_GRID_SIZE = 0.0005
# crystal axis angles in degrees (<001> and <110>)
const CRYSTAL_AXIS_ANGLES = [0, 45]
# SSD adaptive-mesh refinement thresholds (default; can be overridden via metadata settings file)
const DEFAULT_REFINEMENT_LIMITS = [0.2, 0.1, 0.05, 0.02]
# nr of pixels for padding around the map to avoid grid edge effects (default; can be overridden via metadata settings file)
const DEFAULT_PADDING = 3

# Imports for main script
using LegendDataManagement
using LegendHDF5IO
using ArgParse
using PropDicts
using Printf
using Unitful
using LegendSimflow
using Pkg
using SolidStateDetectors
"""
    main()

Generate HPGe drift time maps for specified detector and save to LH5 file.
"""
function main()
    T = Float32
    ssd_version = pkgversion(SolidStateDetectors)
    ldm_version = pkgversion(LegendDataManagement)

    @info "using SSD: $ssd_version"
    @info "using LDM: $ldm_version"

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
        "--ssd-settings"
        help = "Path to ssd settings YAML file (optional; built-in defaults used if absent or missing)"
        default = nothing
    end

    parsed_args = parse_args(s)

    det = parsed_args["detector"]
    meta_path = parsed_args["metadata"]
    opv = parsed_args["opv"]
    output_file = parsed_args["output-file"]
    info_file = parsed_args["info-file"]

    isfile(output_file) && error("Output file already exists")

    raw_opv = parsed_args["opv"]
    opv_val = isnothing(raw_opv) ? nothing : parse(Float32, raw_opv)

    meta, xtal, opv_val = load_detector_metadata(meta_path, det, opv_val)


    # Load optional simulation settings, falling back to built-in defaults.
    # The settings file path is passed via --ssd-settings and applies globally
    # to all detectors and voltages.
    ssd_settings = parsed_args["ssd-settings"]
    sim_cfg = (!isnothing(ssd_settings) && isfile(ssd_settings)) ? readprops(ssd_settings) : PropDict()
    grid_size = get(sim_cfg, :grid_size_in_mm, DEFAULT_GRID_SIZE * 1000) / 1000
    ref_limits = get(sim_cfg, :ssd_refinement_limits, DEFAULT_REFINEMENT_LIMITS)
    padding = get(sim_cfg, :padding, DEFAULT_PADDING)

    sim, info = setup_hpge_simulation(meta_path, meta, xtal, opv_val, T, ref_limits)

    # Compute drift time maps for each crystal axis angle
    output = nothing
    for angle in CRYSTAL_AXIS_ANGLES
        result = compute_drift_time_map(sim, meta, T, angle, grid_size, padding)

        key = Symbol("drift_time_$(lpad(string(angle), 3, '0'))_deg")
        if output === nothing
            output = Dict{Symbol,Any}(pairs(result))
        else
            output[key] = result[key]
        end
    end

    @info "Saving to disk..."
    lh5open(output_file, "cw") do f
        return f[det] = (; output...)
    end

    # the SSD-modeling provenance scalars (impurity scaling, measured and
    # simulated depletion voltages) are stored separately as metadata, not as
    # LH5 scalars; they are aggregated into a `detinfo` file downstream
    if !isnothing(info_file)
        @info "Saving SSD-modeling info to $info_file..."
        writeprops(info_file, PropDict(info))
    end
end


main()
