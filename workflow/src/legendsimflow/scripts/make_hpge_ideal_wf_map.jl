# Script to build ideal waveform map of a detector using ssd
# Map construction based on https://github.com/legend-exp/legend-simflow/blob/main/workflow/src/legendsimflow/scripts/make_hpge_drift_time_maps.jl


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

## Constants
GRID_SIZE = 0.005 # in m
CRYSTAL_AXIS_ANGLES = [0, 45] # in deg <001> <110>
SIM_ENERGY = 2039u"keV"
MAX_NSTEPS = 3000
WAVEFORM_LENGTH = 3000
TIME_STEP = 1u"ns" # sampling rate sim waveform -- T(1)u"ns" ?


## Waveform map construction
function compute_waveform_map_for_angle(
    sim,
    meta,
    T,
    angle_deg,
    only_holes,
    handle_nplus
)
    @info "Computing waveform map at angle $angle_deg deg..."

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
            push!(spawn_positions, CartesianPoint{T}(x * cos(angle_rad), x * sin(angle_rad), z))
            push!(idx_spawn_positions, CartesianIndex(i, k))
        end
    end

    if (!handle_nplus)
        in_idx = findall(
            x -> in(x, sim.detector) && (!in(x, sim.detector.contacts)),
            spawn_positions
        )
    else
        in_idx = findall(x -> in(x, sim.detector), spawn_positions)
    end

    n = length(in_idx)
    wf_signals_threaded = Vector{Vector{Float32}}(undef, n)


    function get_positions(idx, handle_nplus_local = false, verbose = true)
        pos_candidate = spawn_positions[idx]
        if (!handle_nplus_local || (!in(pos_candidate, sim.detector.contacts)))
            return pos_candidate
        else
            min_dist = Inf
            pos_tmp = nothing
            for pos in spawn_positions
                if (!in(pos, sim.detector.contacts) && (in(pos, sim.detector)))
                    dist = norm(pos - pos_candidate)
                    if (dist < min_dist)
                        pos_tmp = pos
                        min_dist = dist
                    end
                end
            end
        end
        return pos_tmp
    end

    @info "Simulating grid (r, z) at angle $(angle_deg)°..."

    @threads for i in 1:n
        if (i % 1000 == 0)
            x = round(100 * i / n)
            println("...simulating $i out of $n ($x %)")
        end

        p = get_positions(in_idx[i], handle_nplus)
        e = SSD.Event([p], [SIM_ENERGY])
        simulate!(e, sim, Δt = TIME_STEP, max_nsteps = MAX_NSTEPS, verbose = false)

        if (only_holes)
            wf = get_electron_and_hole_contribution(e, sim, 1).hole_contribution
        else
            wf = e.waveforms[1]
        end

        wf_signals_threaded[i] = ustrip(wf.signal)
    end

    ## Minimal waveform post-processing
    wf_padded = fill(NaN32, WAVEFORM_LENGTH, length(z_axis), length(x_axis)) # Initialize with NaN. Pixels outside the detector (not in in_idx) will remain NaN.

    for (i, idx) in enumerate(idx_spawn_positions[in_idx])
        raw_signal = wf_signals_threaded[i]

        # Normalize charge waveforms so that their amplitude (energy) = 1
        max_val = maximum(raw_signal)
        if max_val > 0
            signal = raw_signal ./ max_val
        else
            signal = raw_signal
        end

        # Determine actual length to copy from sim (clip if longer than fixed length)
        len = min(length(signal), WAVEFORM_LENGTH)

        # Copy the samples containing the simulated signal
        wf_padded[1:len, idx[2], idx[1]] = signal[1:len]

        # Pad remaining samples with the last value of the charge waveform
        if len < WAVEFORM_LENGTH && len > 0
            last_value = signal[len]
            wf_padded[(len + 1):end, idx[2], idx[1]] .= last_value
        end
    end

    ang_str = lpad(string(angle_deg), 3, '0')

    # Prepare map for writing
    output_fields = Dict{Symbol,Any}(
        :r => collect(x_axis) * u"m",
        :z => collect(z_axis) * u"m",
        :dt => TIME_STEP,
        Symbol("waveform_$(ang_str)_deg") => wf_padded * u"1" # u"1" to indicate dimensionless units (q/q)
    )

    return (; output_fields...)
end


function main()
    SSD = SolidStateDetectors
    T = Float32

    s = ArgParseSettings()

    @add_arg_table s begin
        "--detector"
        help = "HPGe detector name"
        required = true
        "--metadata"
        help = "Path to legend-metadata"
        required = true
        "--output-file"
        help = "Path to output LH5 file"
        required = true
        "--opv"
        help = "detector operational voltage in V (defaults to metadata value)"
        #required = true # Make sure the script is using this!!
        "--use-sqrt" # Outdated??
        help = "Use square-root temperature model"
        action = :store_true
        "--use-sqrt-new" # Outdated??
        help = "Use square-root temperature model with 2016 inputs"
        action = :store_true
        "--only-holes"
        help = "Use only the hole contribution"
        action = :store_true
        "--use-bulk-drift-time"
        help = "Use the drift time for the nearest bulk point for surface events"
        action = :store_true
        "--use-corrections" # Important!! Must be set to true!!
        help = "Use the impurity corrections"
        action = :store_true
    end

    parsed_args = parse_args(s)

    # Extract arguments
    det = parsed_args["detector"]

    meta_path = parsed_args["metadata"]

    output_file = parsed_args["output-file"]
    isfile(output_file) && error("output file already exists")

    # Get the opv - could be a String or Nothing
    raw_opv = parsed_args["opv"]
    if isnothing(raw_opv)
        # 2. If nothing, load the metadata first to find the default
        meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")
        opv_val = Float32(meta.characterization.l200_site.recommended_voltage_in_V)
        @info "No OPV provided. Using metadata default: $opv_val V"
    else
        # If opv provided as keyword argument, parse it
        opv_val = parse(Float32, raw_opv)
        @info "Using user-provided OPV: $opv_val V"
    end

    handle_nplus = parsed_args["use-bulk-drift-time"]
    use_sqrt = parsed_args["use-sqrt"] # Remove?
    only_holes = parsed_args["only-holes"]
    use_corrections = parsed_args["use-corrections"]

    # Load metadata
    meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")
    meta.characterization.l200_site.recommended_voltage_in_V = opv_val

    if (!handle_nplus)
        meta.production.characterization.manufacturer.dl_thickness_in_mm = 0.1
        meta.characterization.combined_0vbb_analysis.fccd_in_mm.value = 0.1
    end

    ids = Dict("bege" => "B", "coax" => "C", "ppc" => "P", "icpc" => "V")
    crystal = ids[meta.type] * @sprintf("%02d", meta.production.order) * meta.production.crystal
    xtal = readprops("$meta_path/hardware/detectors/germanium/crystals/$crystal.yaml")

    if (!use_corrections && (:corrections in keys(xtal.impurity_curve)))
        xtal.impurity_curve.corrections.scale = 1
    end

    # Run ssd simulation
    sim = Simulation{T}(LegendData, meta, xtal,
        HPGeEnvironment("LAr", 87u"K"),
        operational_voltage = opv_val * u"V" #,
        #n_thickness = 0.7u"mm" #????????????????????????????????????????????????????
    )


    sim.detector = SolidStateDetector(
        sim.detector,
        contact_id = 2,
        contact_potential = opv_val
    )


    # Physics Models ## Remove??
    if (use_sqrt)
        charge_drift_model = ADLChargeDriftModel(
            joinpath(
                SolidStateDetectors.get_path_to_example_config_files(),
                "ADLChargeDriftModel",
                "drift_velocity_config_squareroot.yaml"
            )
        )
        sim.detector = SolidStateDetector(sim.detector, charge_drift_model)
    end

    # Calculate potential
    @info "Calculating potential at $(opv_val)V..."
    calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1, 0.05, 0.01], depletion_handling = true)

    @info "Calculating electric field..."
    calculate_electric_field!(sim)

    @info "Calculating depletion voltage..."
    dep = nothing
    try
        dep = estimate_depletion_voltage(sim)
    catch e
        error("Detector is not depleted!")
    end
    @info "Simulated depletion is $dep"


    dep_meas = meta[:characterization][:l200_site][:depletion_voltage_in_V] * u"V"
    @info "Depletion measured during characterization is $dep_meas"

    if abs(dep_meas - dep) > 100 * u"V"
        error("difference between measured and simulated depletion is larger than 100 V!")
    end


    @info "Calculating weighting potential..."
    calculate_weighting_potential!(
        sim,
        sim.detector.contacts[1].id,
        refinement_limits = [0.2, 0.1, 0.05, 0.01],
        verbose = false
    )

    # Compute and save drift-time maps for different angles
    output = Dict{Symbol,Any}()
    for a in CRYSTAL_AXIS_ANGLES
        eff_angle = use_sqrt ? (a - 45) : a
        out = compute_waveform_map_for_angle(sim, meta, T, eff_angle, only_holes, handle_nplus)

        if isempty(output)
            output[:r] = out.r
            output[:z] = out.z
            output[:dt] = out.dt
        end

        key = Symbol("waveform_$(lpad(string(a), 3, '0'))_deg")
        output[key] = out[key]
    end

    @info "Saving to disk..."
    lh5open(output_file, "cw") do f
        return f[det] = (; output...)
    end
end

main()
