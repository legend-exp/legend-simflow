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

# Include helper functions from the library
include(joinpath(@__DIR__, "libjl", "drift_time_helpers.jl"))


"""
    main()

Generate HPGe drift time maps for specified detector and save to LH5 file.
"""
function main()
    SSD = SolidStateDetectors
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

    @info "Operating $det at $opv V"

    meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")
    meta.characterization.l200_site.recommended_voltage_in_V = parse(Float32, opv)

    ids = Dict("bege" => "B", "coax" => "C", "ppc" => "P", "icpc" => "V")
    crystal = ids[meta.type] * @sprintf("%02d", meta.production.order) * meta.production.crystal
    xtal = readprops("$meta_path/hardware/detectors/germanium/crystals/$crystal.yaml")

    # Create simulation
    sim = Simulation{T}(LegendData, meta, xtal)

    charge_drift_model = ADLChargeDriftModel(
        "$meta_path/simprod/config/pars/geds/ssd/adl-2016-temp-model.yaml"
    )
    sim.detector = SolidStateDetector(sim.detector, charge_drift_model)
    sim.detector = SolidStateDetector(
        sim.detector,
        contact_id = 2,
        contact_potential = meta.characterization.l200_site.recommended_voltage_in_V
    )

    @info "Calculating electric potential..."
    calculate_electric_potential!(
        sim,
        refinement_limits = REFINEMENT_LIMITS,
        depletion_handling = true
    )

    @info "Calculating electric field..."
    calculate_electric_field!(sim)

    dep = nothing
    try
        dep = estimate_depletion_voltage(sim)
    catch
        error("Detector is not depleted!")
    end
    @info "Simulated depletion is $dep"

    dep_meas = meta[:characterization][:l200_site][:depletion_voltage_in_V] * u"V"
    @info "Depletion measured during characterization is $dep_meas"

    if abs(dep_meas - dep) > 100 * u"V"
        error("Difference between measured and simulated depletion is larger than 100 V!")
    end

    @info "Calculating weighting potential..."
    calculate_weighting_potential!(
        sim,
        sim.detector.contacts[1].id,
        refinement_limits = REFINEMENT_LIMITS,
        verbose = false
    )

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
