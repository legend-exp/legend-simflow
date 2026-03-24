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

# grid spacing in meters
const GRID_SIZE = 0.0005
# crystal axis angles in degrees (<001> and <110>)
const CRYSTAL_AXIS_ANGLES = [0, 45]
# SSD adaptive-mesh refinement thresholds as fractions of the crystal radius
# matches current SSD behaviour
const REFINEMENT_LIMITS = [0.2, 0.1, 0.05, 0.02]
# nr of pixels for padding around the map to avoid grid edge effects
const PADDING = 3

# Imports for main script
using LegendDataManagement
using LegendHDF5IO
using ArgParse
using PropDicts
using Printf
using Unitful
using LegendSimflow

"""
    main()

Generate HPGe drift time maps for specified detector and save to LH5 file.
"""
function main()
    T = Float32

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
        "--opv"
        help = "detector operational voltage in V"
        required = true
    end
    @add_arg_table s begin
        "--output-file"
        help = "Path to output LH5 file"
        required = true
    end

    parsed_args = parse_args(s)

    det = parsed_args["detector"]
    meta_path = parsed_args["metadata"]
    opv = parsed_args["opv"]
    output_file = parsed_args["output-file"]

    isfile(output_file) && error("Output file already exists")

    raw_opv = parsed_args["opv"]
    opv_val = isnothing(raw_opv) ? nothing : parse(Float32, raw_opv)

    meta, xtal, opv_val = load_detector_metadata(meta_path, det, opv_val)

    sim = setup_hpge_simulation(meta_path, meta, xtal, opv_val, T, REFINEMENT_LIMITS)

    # Compute drift time maps for each crystal axis angle
    output = nothing
    for angle in CRYSTAL_AXIS_ANGLES
        result = compute_drift_time_map(sim, meta, T, angle, GRID_SIZE, PADDING)

        key = Symbol("drift_time_$(lpad(string(angle), 3, '0'))_deg")
        if output === nothing
            output = Dict(pairs(result))
        else
            output[key] = result[key]
        end
    end

    @info "Saving to disk..."
    lh5open(output_file, "cw") do f
        return f[det] = (; output...)
    end
end


main()
