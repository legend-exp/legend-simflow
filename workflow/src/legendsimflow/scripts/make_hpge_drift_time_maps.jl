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

using SolidStateDetectors
using LegendDataManagement
using Unitful
using LegendHDF5IO
using Base.Threads
using LinearAlgebra
using ArgParse
using PropDicts
using Printf
using RadiationDetectorDSP


GRID_SIZE = 0.0005 # in m
CRYSTAL_AXIS_ANGLES = [0, 45] # in deg


function compute_drift_time(wf, rise_convergence_criteria, tint)
    collected_charge = wf[argmax(abs.(wf))]

    if collected_charge < 0 # very rare, occurs when electron drift dominates and holes are stuck
        wf *= -1 # to ensure Intersect() works as intended
        collected_charge *= -1
    end
    intersection = tint(wf, rise_convergence_criteria * collected_charge)
    dt_intersection = ceil(Int, intersection.x)
    dt_fallback = length(wf)
    dt_diff = dt_fallback - dt_intersection

    dt = if intersection.multiplicity > 0
        if dt_diff > 2 # dt_intersection is not at the end of the wf
            tint2 = Intersect(mintot = dt_diff) # check intersect again but with min_n_over_thresh (mintot) set to max
            intersection2 = tint2(wf, rise_convergence_criteria * collected_charge)
            if intersection2.multiplicity > 0  # usually monotonic waveforms which converge very slowly
                dt_intersection
            else # usually non-monotonic waveforms
                dt_fallback
            end
        else # dt_intersection is at the end of the wf (ie. both drift length and wf convergence methods agree)
            dt_intersection
        end
    else # no intersection with minimal conditions (mintot=0) found
        dt_fallback
    end
    return dt
end


function compute_drift_map_for_angle(
    sim,
    meta,
    T,
    angle_deg
)
    @info "Computing drift time map at angle $angle_deg deg..."

    SSD = SolidStateDetectors
    angle_rad = deg2rad(angle_deg)

    function make_axis(T, boundary, GRID_SIZE)
        # define interior domain strictly within (0, boundary)
        offset = 2 * SSD.ConstructiveSolidGeometry.csg_default_tol(T)
        inner_start = 0 + offset
        inner_stop = boundary - offset

        # compute number of intervals in the interior (ensure at least 1)
        n = max(1, round(Int, (inner_stop - inner_start) / GRID_SIZE))

        # recompute step to fit the inner domain evenly
        step = (inner_stop - inner_start) / n

        # create interior axis
        axis = range(inner_start, step = step, length = n + 1)

        # prepend and append slightly out-of-bound points
        extended_axis = [0 - offset, axis..., boundary + offset]

        return extended_axis
    end

    radius = meta.geometry.radius_in_mm / 1000
    height = meta.geometry.height_in_mm / 1000

    x_axis = make_axis(T, radius, GRID_SIZE)
    z_axis = make_axis(T, height, GRID_SIZE)

    spawn_positions = CartesianPoint{T}[]
    idx_spawn_positions = CartesianIndex[]

    for (i, x) in enumerate(x_axis)
        for (k, z) in enumerate(z_axis)
            point = T[x * cos(angle_rad), x * sin(angle_rad), z]
            push!(spawn_positions, CartesianPoint(point))
            push!(idx_spawn_positions, CartesianIndex(i, k))
        end
    end
    in_idx = findall(x -> in(x, sim.detector), spawn_positions)

    # simulate events
    time_step = T(1)u"ns"
    max_nsteps = 10000

    # prepare thread-local storage
    n = length(in_idx)
    dt_threaded = Vector{Int}(undef, n)
    rise_convergence_criteria = 1 - 1e-6
    tint = Intersect(mintot = 0)

    function deriv(wf)
        return argmax(ustrip(wf.signal))
    end

    function get_positions(idx, verbose = true)
        pos_candidate = spawn_positions[idx]

        # if the point isn't inside the contact do nothing
        if ((!in(pos_candidate, sim.detector.contacts)))
            return pos_candidate
        else
            min_dist = Inf
            pos_tmp = nothing

            if (verbose)
                @debug "position $pos_candidate is in the contact searching for a new one"
            end

            for pos in spawn_positions
                if (!in(pos, sim.detector.contacts) && (in(pos, sim.detector)))
                    dist = norm(pos - pos_candidate)

                    if (dist < min_dist)
                        pos_tmp = pos
                        min_dist = dist
                    end
                end
            end
            if (verbose)
                @debug "found position $pos_tmp $min_dist away"
            end
        end
        return pos_tmp
    end

    @info "Simulating energy depositions on grid r=0:$GRID_SIZE:$radius and z=0:$GRID_SIZE:$height at angle $(angle_deg)°..."
    @threads for i in 1:n
        if (i % 1000 == 0)
            x = round(100 * i / n)
            println("...simulating $i out of $n ($x %)")
        end

        p = get_positions(in_idx[i])
        e = SSD.Event([p], [2039u"keV"])
        simulate!(e, sim, Δt = time_step, max_nsteps = max_nsteps, verbose = false)

        wf = get_electron_and_hole_contribution(e, sim, 1).hole_contribution

        dt_threaded[i] =
            compute_drift_time(ustrip(wf.signal), rise_convergence_criteria, tint)
    end
    dt = dt_threaded

    drift_time = fill(NaN, length(x_axis), length(z_axis))
    for (i, idx) in enumerate(idx_spawn_positions[in_idx])
        drift_time[idx] = dt[i]
    end

    ang_str = lpad(string(angle_deg), 3, '0')
    output = (;
        :r => collect(x_axis) * u"m",
        :z => collect(z_axis) * u"m",
        Symbol("drift_time_$(ang_str)_deg") => transpose(drift_time) * u"ns"
    )

    return output
end

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

    isfile(output_file) && error("output file already exists")

    meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")

    meta.characterization.l200_site.recommended_voltage_in_V = parse(Float32, opv)

    ids = Dict("bege" => "B", "coax" => "C", "ppc" => "P", "icpc" => "V")
    crystal =
        ids[meta.type] * @sprintf("%02d", meta.production.order) * meta.production.crystal
    xtal = readprops("$meta_path/hardware/detectors/germanium/crystals/$crystal.yaml")

    # make the simulation
    sim = Simulation{T}(LegendData, meta, xtal)

    charge_drift_model = ADLChargeDriftModel(
        "$meta_path/simprod/config/pars/ssd/adl-2016-temp-model.yaml"
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
        refinement_limits = [0.2, 0.1, 0.05, 0.01],
        depletion_handling = true
    )

    @info "Calculating electric field..."
    calculate_electric_field!(sim)
    dep = nothing
    try
        dep = estimate_depletion_voltage(sim)
    catch e
        @info "detector is not depleted"
    end

    dep_meas = meta[:characterization][:l200_site][:depletion_voltage_in_V]

    @info "Measured depletion is $dep_meas V"
    @info "Simulated depletion is $dep V"
    @info "Calculating weighting potential..."
    calculate_weighting_potential!(
        sim,
        sim.detector.contacts[1].id,
        refinement_limits = [0.2, 0.1, 0.05, 0.01],
        verbose = false
    )

    # Compute and save drift-time maps for different angles
    output = nothing
    for a in CRYSTAL_AXIS_ANGLES
        eff_angle = use_sqrt ? (a - 45) : a
        out = compute_drift_map_for_angle(
            sim,
            meta,
            T,
            eff_angle
        )

        key = Symbol("drift_time_$(lpad(string(a), 3, '0'))_deg")
        if output == nothing
            output = Dict(pairs(out))
        else
            output[key] = out[key]
        end

    end

    @info "Saving to disk..."
    lh5open(output_file, "cw") do f
        return f[det] = (; output...)
    end
end


main()
