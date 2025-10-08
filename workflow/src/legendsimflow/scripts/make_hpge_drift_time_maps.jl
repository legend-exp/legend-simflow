using SolidStateDetectors
using LegendTestData
using LegendDataManagement
using Unitful
using LegendHDF5IO
using Base.Threads
using LinearAlgebra
using ArgParse
using PropDicts
using Printf
using RadiationDetectorDSP


function compute_drift_time(wf, rise_convergence_creteria, tint)
    collected_charge = wf[argmax(abs.(wf))]

    if collected_charge < 0 # very rare, occurs when electron drift dominates and holes are stuck
        wf *= -1 # to ensure Intersect() works as intended
        collected_charge *= -1
    end
    intersection = tint(wf, rise_convergence_creteria * collected_charge)
    dt_intersection = ceil(intersection.x)
    dt_fallback = length(wf)
    dt_diff = dt_fallback - dt_intersection

    dt = if intersection.multiplicity > 0
        if dt_diff > 2 # dt_intersection is not at the end of the wf
            tint2 = Intersect(mintot = dt_diff) # check intersect again but with min_n_over_thresh (mintot) set to max
            intersection2 = tint2(wf, rise_convergence_creteria * collected_charge)
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
    dt
end


function main()
    SSD = SolidStateDetectors
    T = Float32

    s = ArgParseSettings()

    @add_arg_table s "--detector" help = "HPGe detector name" required = true
    @add_arg_table s "--metadata" help = "Path to legend-metadata" required = true
    @add_arg_table s "--opv-file" help = "Path to file with operational voltages" required =
        true
    @add_arg_table s "--output-file" help = "Path to output LH5 file" required = true

    @add_arg_table s "--angle" help = "Rotation angle" default = "0"

    @add_arg_table s begin
        "--use-sqrt"
        help = "Use square-root temperature model"
        action = :store_true
    end
    @add_arg_table s begin
        "--use-sqrt-new"
        help = "Use square-root temperature model with 2016 inputs"
        action = :store_true
    end
    @add_arg_table s begin
        "--get-max-time"
        help = "Get the time of the maximum"
        action = :store_true
    end
    @add_arg_table s begin
        "--only-holes"
        help = "Use only the hole contribution"
        action = :store_true
    end
    @add_arg_table s begin
        "--use-bulk-drift-time"
        help = "Use the drift time for the nearest bulk point for surface events"
        action = :store_true
    end
    @add_arg_table s begin
        "--use-corrections"
        help = "Use the corrections"
        action = :store_true
    end

    get_drift_time = true
    use_triangle = false
    use_sqrt = false
    use_sqrt_new = true
    handle_nplus = true
    only_holes = true
    use_corrections = true

    parsed_args = parse_args(s)

    det = parsed_args["detector"]
    meta_path = parsed_args["metadata"]
    opv_path = parsed_args["opv-file"]
    output_file = parsed_args["output-file"]

    use_sqrt = parsed_args["use-sqrt"]
    handle_nplus = parsed_args["use-bulk-drift-time"]
    use_sqrt_new = parsed_args["use-sqrt-new"]
    get_drift_time = !parsed_args["get-max-time"]
    only_holes = parsed_args["only-holes"]
    angle = parse(T, parsed_args["angle"])
    use_corrections = parsed_args["use-corrections"]

    angle_orig = angle
    if use_sqrt
        angle = angle - 45
    end

    nphi = 36
    if angle == 0.0
        nphi = 2
    end
    angle_rad = deg2rad(angle)

    meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")

    meta.characterization.l200_site.recommended_voltage_in_V =
        readprops("$opv_path")[Symbol(det)].operational_voltage_in_V

    if (!handle_nplus)
        meta.production.characterization.manufacturer.dl_thickness_in_mm = 0.1
        meta.characterization.combined_0vbb_analysis.fccd_in_mm.value = 0.1
    end

    crystal = "V" * @sprintf("%02d", meta.production.order) * meta.production.crystal
    xtal = readprops("$meta_path/hardware/detectors/germanium/crystals/$crystal.yaml")

    if (!use_corrections && (:corrections in keys(xtal.impurity_curve)))
        xtal.impurity_curve.corrections.scale = 1
    end

    # make the simulation
    sim = Simulation{T}(LegendData, meta, xtal)

    if (use_sqrt)
        charge_drift_model = ADLChargeDriftModel(
            joinpath(
                SolidStateDetectors.get_path_to_example_config_files(),
                "ADLChargeDriftModel",
                "drift_velocity_config_squareroot.yaml",
            ),
        )
        sim.detector = SolidStateDetector(sim.detector, charge_drift_model)
    end

    if (use_sqrt_new)
        charge_drift_model = ADLChargeDriftModel(
            "$meta_path/simprod/config/pars/ssd/adl-2016-temp-model.yaml",
        )
        sim.detector = SolidStateDetector(sim.detector, charge_drift_model)
    end

    sim.detector = SolidStateDetector(
        sim.detector,
        contact_id = 2,
        contact_potential = meta.characterization.l200_site.recommended_voltage_in_V,
    )

    @info "Calculating electric potential..."
    calculate_electric_potential!(
        sim,
        refinement_limits = [0.2, 0.1, 0.05, 0.01],
        depletion_handling = true,
    )

    @info "Calculating electric field..."
    calculate_electric_field!(sim, n_points_in_φ = nphi)
    dep = nothing
    try
        dep = estimate_depletion_voltage(sim)
    catch e
        @info "detector is not depleted"
    end

    dep_meas = meta[:characterization][:l200_site][:depletion_voltage_in_V]

    @info "Measured depletion is $dep_meas V"
    @info "Simulationed depletion = $dep V"
    @info "Calculating weighting potential..."
    calculate_weighting_potential!(
        sim,
        sim.detector.contacts[1].id,
        refinement_limits = [0.2, 0.1, 0.05, 0.01],
        verbose = false,
        n_points_in_φ = nphi,
    )

    function make_axis(T, boundary, gridsize)
        # define interior domain strictly within (0, boundary)
        offset = 2 * SSD.ConstructiveSolidGeometry.csg_default_tol(T)
        inner_start = 0 + offset
        inner_stop = boundary - offset

        # compute number of intervals in the interior
        n = round(Int, (inner_stop - inner_start) / gridsize)

        # recompute step to fit the inner domain evenly
        step = (inner_stop - inner_start) / n

        # create interior axis
        axis = range(inner_start, step = step, length = n + 1)

        # prepend and append slightly out-of-bound points
        extended_axis = [0 - offset, axis..., boundary + offset]

        return extended_axis
    end

    gridsize = 0.0005 # in m
    radius = meta.geometry.radius_in_mm / 1000
    height = meta.geometry.height_in_mm / 1000

    x_axis = make_axis(T, radius, gridsize)
    z_axis = make_axis(T, height, gridsize)

    spawn_positions = CartesianPoint{T}[]
    idx_spawn_positions = CartesianIndex[]

    for (i, x) in enumerate(x_axis)
        for (k, z) in enumerate(z_axis)
            point = T[x*cos(angle_rad), x*sin(angle_rad), z]
            push!(spawn_positions, CartesianPoint(point))
            push!(idx_spawn_positions, CartesianIndex(i, k))
        end
    end
    if (!handle_nplus)
        in_idx = findall(
            x -> in(x, sim.detector) && (!in(x, sim.detector.contacts)),
            spawn_positions,
        )
    else
        in_idx = findall(x -> in(x, sim.detector), spawn_positions)

    end

    # simulate events

    time_step = T(1)u"ns"
    max_nsteps = 10000

    # prepare thread-local storage
    n = length(in_idx)
    dt_threaded = Vector{Int}(undef, n)
    rise_convergence_creteria = 1-1e-6
    tint = Intersect(mintot = 0)

    function deriv(wf)
        if (!use_triangle)
            return argmax(diff(ustrip(wf.signal)))
        else
            trap = TrapezoidalChargeFilter(50u"ns", 5u"ns", 50u"ns");

            wf_shift = RDWaveform(wf.time, wf.signal .- wf.signal[begin])
            wf_ext = add_baseline_and_extend_tail(wf_shift, 100, 4000)
            wf_no_units = RDWaveform(wf_ext.time, ustrip.(wf_ext.signal))
            curr_wf = trap(wf_no_units)

            idx_A = argmax(curr_wf.signal)
            return currwf.time[idx_A]
        end
    end

    function get_positions(idx, handle_nplus = false, verbose = true)
        pos_candidate = spawn_positions[idx]
        # if we dont shift points inside or if the point isn't inside the contact do nothing
        if (!handle_nplus | (!in(pos_candidate, sim.detector.contacts)))
            return pos_candidate
        else
            min_dist = Inf
            pos_tmp = nothing

            if (verbose)
                @debug "position $pos_candidate is in the contact searching for a new one"
            end

            for pos in spawn_positions

                if (!in(pos, sim.detector.contacts) & (in(pos, sim.detector)))
                    dist = norm(pos - pos_candidate)

                    if (dist<min_dist)
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

    @info "Simulating energy depositions on grid r=0:$gridsize:$radius and z=0:$gridsize:$height..."
    @threads for i = 1:n
        if (i % 1000 == 0)
            x = round(100 * i / n)
            println("...simulating $i out of $n ($x %)")
        end

        p = get_positions(in_idx[i], handle_nplus)
        e = SSD.Event([p], [2039u"keV"])
        simulate!(e, sim, Δt = time_step, max_nsteps = max_nsteps, verbose = false)

        if (only_holes)
            wf = get_electron_and_hole_contribution(e, sim, 1).hole_contribution
        else
            wf = e.waveforms[1]
        end

        if (get_drift_time)
            dt_threaded[i] =
                compute_drift_time(ustrip(wf.signal), rise_convergence_creteria, tint)
        else
            dt_threaded[i] = deriv(wf)

        end
    end
    dt = dt_threaded

    drift_time = fill(NaN, length(x_axis), length(z_axis))
    for (i, idx) in enumerate(idx_spawn_positions[in_idx])
        drift_time[idx] = dt[i]
    end

    output = (
        r = collect(x_axis) * u"m",
        z = collect(z_axis) * u"m",
        drift_time = transpose(drift_time) * u"ns",
    )

    @info "Saving to disk..."

    lh5open(output_file, "w") do f
        f[det] = output
    end
end


main()
