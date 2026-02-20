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

using SolidStateDetectors
using LegendDataManagement
using Unitful
using LegendHDF5IO
using Base.Threads
using ArgParse
using PropDicts
using Printf
using RadiationDetectorDSP

include(joinpath(@__DIR__, "libjl", "drift_time_helpers.jl"))

## Constants
GRID_SIZE = 0.0005 # cylindrical (r, z) grid spacing in meters
CRYSTAL_AXIS_ANGLES = [0, 45] # crystal axis orientation angles in degrees (<001>, <110>)
SIM_ENERGY = 2039u"keV" # mono-energetic deposition used to generate ideal waveforms
MAX_NSTEPS = 3000 # maximum number of integration steps per drift trajectory (dimensionless)
WAVEFORM_LENGTH = 3000 # number of time samples per simulated waveform (dimensionless)
TIME_STEP = 1u"ns" # time step between waveform samples


## Waveform map construction
"""
    compute_waveform_map_for_angle(
        sim::SolidStateDetectors.Simulation,
        meta,
        T::Type{<:AbstractFloat},
        angle_deg::Real,
        only_holes::Bool,
        handle_nplus::Bool
    ) -> NamedTuple

Compute a 2D waveform map (r, z grid) at a specified azimuthal angle.

Simulates charge carrier drift for a grid of positions in cylindrical coordinates,
generating waveforms that represent the detector response. The grid is constructed
in the (r, z) plane and then rotated to the specified angle.

# Arguments
- `sim::SolidStateDetectors.Simulation`: The configured detector simulation object
- `meta`: Metadata object containing detector geometry (radius_in_mm, height_in_mm)
- `T::Type{<:AbstractFloat}`: Floating-point precision type (typically Float32)
- `angle_deg::Real`: Azimuthal angle in degrees for the r-z plane rotation
- `only_holes::Bool`: If true, extract only hole contribution; if false, use full waveform
- `handle_nplus::Bool`: If true, handle positions at n+ contact by finding nearest valid point

# Returns
- `NamedTuple`: Contains the following fields:
  - `r`: Radial coordinates (m) as a vector with units
  - `z`: Axial coordinates (m) as a vector with units
  - `dt`: Time step (sampling rate) with units
  - `waveform_XXX_deg`: 3D array [time_samples, n_z, n_r] of normalized waveforms (dimensionless)
    where XXX is the zero-padded angle (e.g., waveform_045_deg for 45°)

# Notes
- Grid resolution is controlled by the global constant `GRID_SIZE`
- Waveform length is controlled by `WAVEFORM_LENGTH`
- Positions outside the detector remain as NaN in the output array
- All waveforms are normalized so max(waveform) = 1.0
- Grid construction and nearest-bulk-point lookup use helper functions from `drift_time_helpers.jl`
"""
function compute_waveform_map_for_angle(
    sim::SolidStateDetectors.Simulation,
    meta,
    T::Type{<:AbstractFloat},
    angle_deg::Real,
    only_holes::Bool,
    handle_nplus::Bool
)
    @info "Computing waveform map at angle $angle_deg deg..."

    SSD = SolidStateDetectors
    angle_rad = deg2rad(angle_deg)

    radius = meta.geometry.radius_in_mm / 1000
    height = meta.geometry.height_in_mm / 1000

    # Use functions from libjl/drift_time_helpers.jl
    x_axis = build_simulation_grid_axis(T, radius, GRID_SIZE)
    z_axis = build_simulation_grid_axis(T, height, GRID_SIZE)

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

    @info "Simulating grid (r, z) at angle $(angle_deg)°..."

    @threads for i in 1:n
        if (i % 1000 == 0)
            x = round(100 * i / n)
            @info "...simulating $i out of $n ($x %)"
        end

        # Use functions from libjl/drift_time_helpers.jl
        p = find_valid_spawn_position(in_idx[i], spawn_positions, sim.detector; verbose = false)
        if p === nothing
            @warn "find_valid_spawn_position did not return a valid spawn position for index $(in_idx[i]); skipping."
            continue
        end


        e = SSD.Event([p], [SIM_ENERGY])
        simulate!(e, sim, Δt = TIME_STEP, max_nsteps = MAX_NSTEPS, verbose = false)

        if (only_holes)
            wf = get_electron_and_hole_contribution(e, sim, 1).hole_contribution
        else
            wf = e.waveforms[1]
        end

        wf_signals_threaded[i] = ustrip(wf.signal)
    end

    # Minimal waveform post-processing
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
        "--use-corrections"
        help = "Apply impurity profile corrections from crystal characterization (default: enabled)"
        action = :store_true
        default = true
    end
    @add_arg_table s begin
        "--opv"
        help = "detector operational voltage in V (defaults to metadata value)"
    end
    @add_arg_table s begin
        "--use-sqrt" # REVIEW: Outdated??
        help = "Use square-root temperature model"
        action = :store_true
    end
    @add_arg_table s begin
        "--use-sqrt-new" # REVIEW: Outdated??
        help = "Use square-root temperature model with 2016 inputs"
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

    parsed_args = parse_args(s)

    # Extract arguments
    det = parsed_args["detector"]

    meta_path = parsed_args["metadata"]

    output_file = parsed_args["output-file"]

    # Let the file be overwritten if it exists (matches Python behavior) # REVIEW: Allow overwriting?
    if isfile(output_file)
        @info "Output file exists and will be overwritten: $output_file"
    end


    # Load metadata
    meta = readprops("$meta_path/hardware/detectors/germanium/diodes/$det.yaml")

    # Get the opv - could be a String or Nothing
    raw_opv = parsed_args["opv"]
    if isnothing(raw_opv)
        # If nothing, use the metadata value
        opv_val = Float32(meta.characterization.l200_site.recommended_voltage_in_V)
        @info "No OPV provided. Using metadata default: $opv_val V"
    else
        # If opv provided as keyword argument, parse it
        opv_val = parse(Float32, raw_opv)
        @info "Using user-provided OPV: $opv_val V"
    end
    meta.characterization.l200_site.recommended_voltage_in_V = opv_val

    handle_nplus = parsed_args["use-bulk-drift-time"]
    use_sqrt = parsed_args["use-sqrt"]
    only_holes = parsed_args["only-holes"]
    use_corrections = parsed_args["use-corrections"]
    if !use_corrections
        @warn "Impurity corrections disabled! This may produce unrealistic results."
    end




    if (!handle_nplus)
        meta.production.characterization.manufacturer.dl_thickness_in_mm = 0.1
        meta.characterization.combined_0vbb_analysis.fccd_in_mm.value = 0.1
    end

    ids = Dict("bege" => "B", "coax" => "C", "ppc" => "P", "icpc" => "V")
    crystal = ids[meta.type] * @sprintf("%02d", meta.production.order) * meta.production.crystal
    xtal = readprops("$meta_path/hardware/detectors/germanium/crystals/$crystal.yaml")

    if !use_corrections
        if :corrections in keys(xtal.impurity_curve)
            xtal.impurity_curve.corrections.scale = 1
        else
            @warn "Impurity corrections were requested to be disabled, but no corrections are defined for crystal $crystal."
        end
    end

    # Run ssd simulation
    sim = Simulation{T}(LegendData, meta, xtal,
        HPGeEnvironment("LAr", 87u"K"),
        operational_voltage = opv_val * u"V"
    )


    sim.detector = SolidStateDetector(
        sim.detector,
        contact_id = 2,
        contact_potential = opv_val
    )


    # Physics Models # REVIEW: Remove??
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

    # Compute and save ideal waveforms map for different angles
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
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        @info "Creating output directory: $output_dir"
        mkpath(output_dir)
    end

    lh5open(output_file, "w") do f
        return f[det] = (; output...)
    end
end

main()
