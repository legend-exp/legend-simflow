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

# Imports needed by this library (not needed in main script)
using SolidStateDetectors
using Unitful
using Base.Threads
using LinearAlgebra
using RadiationDetectorDSP

# 8-connectivity neighbor offsets for map extension
const NEIGHBOR_OFFSETS_8CONN = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


"""
    extract_drift_time_from_waveform(wf, convergence_threshold, intersect_op)

Extract the drift time from a charge waveform using intersection-based analysis.

# Arguments
- `wf`: Waveform signal array (unitless)
- `convergence_threshold`: Fraction of collected charge to define convergence (e.g., 1 - 1e-6)
- `intersect_op`: An `Intersect` operator from RadiationDetectorDSP

# Returns
- `Int`: Drift time in samples
"""
function extract_drift_time_from_waveform(wf, convergence_threshold, intersect_op)
    collected_charge = wf[argmax(abs.(wf))]

    # Handle rare case where electron drift dominates and holes are stuck
    if collected_charge < 0
        wf = wf .* -1  # Create negated copy for Intersect compatibility
        collected_charge *= -1
    end

    threshold_level = convergence_threshold * collected_charge
    intersection = intersect_op(wf, threshold_level)
    dt_intersection = ceil(Int, intersection.x)
    dt_fallback = length(wf)
    dt_diff = dt_fallback - dt_intersection

    dt = if intersection.multiplicity > 0
        if dt_diff > 2
            # Check intersection again with stricter mintot
            intersect_recheck = Intersect(mintot = dt_diff)
            intersection2 = intersect_recheck(wf, threshold_level)
            if intersection2.multiplicity > 0
                # Monotonic waveforms converging slowly
                dt_intersection
            else
                # Non-monotonic waveforms
                dt_fallback
            end
        else
            # Intersection at waveform end (drift length and convergence agree)
            dt_intersection
        end
    else
        # No intersection found
        dt_fallback
    end

    return dt
end


"""
    extend_drift_time_map(drift_map, row_axis, col_axis; layers=1)

Extend a drift time map by adding pixel layers around the valid (non-NaN) region.
New pixel values are set to the average of neighboring non-NaN pixels (8-connectivity).
The grid is internally enlarged by `layers` on applicable sides to ensure the extended map is not clipped.
Note: The column axis (typically radius r) is not extended into negative values.

# Arguments
- `drift_map::AbstractMatrix`: 2D matrix where NaN indicates invalid regions (rows × cols)
- `row_axis::AbstractVector`: Axis values corresponding to rows (e.g., z height)
- `col_axis::AbstractVector`: Axis values corresponding to columns (e.g., r radius, must be non-negative)
- `layers::Int`: Number of pixel layers to add (default: 1)

# Returns
- `NamedTuple`: Contains extended `:drift_map`, `:row_axis`, and `:col_axis`
"""
function extend_drift_time_map(
    drift_map::AbstractMatrix,
    row_axis::AbstractVector,
    col_axis::AbstractVector;
    layers::Int = 1
)
    orig_nrows, orig_ncols = size(drift_map)

    # Determine extension on each side
    # Row axis (z): extend on both sides
    # Col axis (r): only extend on high side (don't go negative)
    row_layers_low = layers
    row_layers_high = layers
    col_layers_low = 0  # Don't extend r into negative values
    col_layers_high = layers

    # Create enlarged grid
    new_nrows = orig_nrows + row_layers_low + row_layers_high
    new_ncols = orig_ncols + col_layers_low + col_layers_high
    work_map = fill(NaN, new_nrows, new_ncols)

    # Copy original data into the grid (offset by the low-side layers)
    work_map[(row_layers_low + 1):(row_layers_low + orig_nrows), (col_layers_low + 1):(col_layers_low + orig_ncols)] =
        drift_map

    # Extend the non-NaN region layer by layer
    for _ in 1:layers
        cur_nrows, cur_ncols = size(work_map)
        boundary_pixels = Tuple{Int,Int}[]
        boundary_values = Float64[]

        for row in 1:cur_nrows, col in 1:cur_ncols
            if isnan(work_map[row, col])
                neighbor_vals = Float64[]
                for (dr, dc) in NEIGHBOR_OFFSETS_8CONN
                    nr, nc = row + dr, col + dc
                    if 1 <= nr <= cur_nrows && 1 <= nc <= cur_ncols
                        val = work_map[nr, nc]
                        !isnan(val) && push!(neighbor_vals, val)
                    end
                end

                if !isempty(neighbor_vals)
                    push!(boundary_pixels, (row, col))
                    push!(boundary_values, sum(neighbor_vals) / length(neighbor_vals))
                end
            end
        end

        # Apply new values after scanning to avoid affecting current layer
        for (idx, (row, col)) in enumerate(boundary_pixels)
            work_map[row, col] = boundary_values[idx]
        end
    end

    # Extend axes to match the enlarged grid
    row_vals = ustrip.(row_axis)
    col_vals = ustrip.(col_axis)

    row_step = maximum(abs.(diff(row_vals)))
    col_step = maximum(abs.(diff(col_vals)))

    # Extend row axis on both sides
    extended_row_vals = vcat(
        [row_vals[1] - row_step * i for i in row_layers_low:-1:1],
        row_vals,
        [row_vals[end] + row_step * i for i in 1:row_layers_high]
    )

    # Extend col axis only on high side (no negative r values)
    extended_col_vals = vcat(
        col_vals,
        [col_vals[end] + col_step * i for i in 1:col_layers_high]
    )

    return (;
        drift_map = work_map,
        row_axis = extended_row_vals * unit(eltype(row_axis)),
        col_axis = extended_col_vals * unit(eltype(col_axis))
    )
end


"""
    build_simulation_grid_axis(T, boundary, grid_step)

Build a grid axis for detector simulation with extended boundary points.

# Arguments
- `T`: Numeric type (e.g., Float32)
- `boundary`: Maximum value of the axis (e.g., detector radius or height in meters)
- `grid_step`: Grid spacing in meters

# Returns
- `Vector`: Grid axis including slightly out-of-bound points at both ends
"""
function build_simulation_grid_axis(T, boundary, grid_step)
    SSD = SolidStateDetectors
    offset = 2 * SSD.ConstructiveSolidGeometry.csg_default_tol(T)
    inner_start = 0 + offset  # Interior domain starts strictly within (0, boundary)
    inner_stop = boundary - offset

    n_intervals = max(1, round(Int, (inner_stop - inner_start) / grid_step))
    step = (inner_stop - inner_start) / n_intervals
    interior_axis = range(inner_start, step = step, length = n_intervals + 1)

    return [-offset, interior_axis..., boundary + offset]
end


"""
    find_valid_spawn_position(candidate_idx, spawn_positions, detector; verbose=true)

Find a valid spawn position outside detector contacts. If the candidate position
is inside a contact, find the nearest valid position inside the detector.

# Arguments
- `candidate_idx`: Index into `spawn_positions` array
- `spawn_positions`: Array of CartesianPoint positions
- `detector`: SolidStateDetector object

# Returns
- `CartesianPoint`: Valid position for event simulation
"""
function find_valid_spawn_position(candidate_idx, spawn_positions, detector; verbose = true)
    pos_candidate = spawn_positions[candidate_idx]

    if !in(pos_candidate, detector.contacts)
        return pos_candidate
    end

    verbose && @debug "Position $pos_candidate is in contact, searching for alternative"

    min_dist = Inf
    pos_result = nothing

    for pos in spawn_positions
        if !in(pos, detector.contacts) && in(pos, detector)
            dist = norm(pos - pos_candidate)
            if dist < min_dist
                pos_result = pos
                min_dist = dist
            end
        end
    end

    verbose && @debug "Found position $pos_result at distance $min_dist"

    return pos_result
end


"""
    compute_drift_time_map(sim, meta, T, angle_deg, grid_step)

Compute a drift time map for an HPGe detector at a specific crystal axis angle.

# Arguments
- `sim`: SolidStateDetectors Simulation object
- `meta`: Detector metadata (PropDict)
- `T`: Numeric type (e.g., Float32)
- `angle_deg`: Crystal axis angle in degrees
- `grid_step`: Grid spacing in meters

# Returns
- `NamedTuple`: Contains `:r`, `:z` axes and `:drift_time_XXX_deg` matrix
"""
function compute_drift_time_map(sim, meta, T, angle_deg, grid_step)
    @info "Computing drift time map at angle $angle_deg deg..."

    SSD = SolidStateDetectors
    angle_rad = deg2rad(angle_deg)

    radius = meta.geometry.radius_in_mm / 1000
    height = meta.geometry.height_in_mm / 1000

    r_axis = build_simulation_grid_axis(T, radius, grid_step)
    z_axis = build_simulation_grid_axis(T, height, grid_step)

    # Build spawn positions grid
    spawn_positions = CartesianPoint{T}[]
    idx_spawn_positions = CartesianIndex[]

    for (i, r) in enumerate(r_axis), (k, z) in enumerate(z_axis)
        point = T[r * cos(angle_rad), r * sin(angle_rad), z]
        push!(spawn_positions, CartesianPoint(point))
        push!(idx_spawn_positions, CartesianIndex(i, k))
    end

    inside_detector_idx = findall(p -> in(p, sim.detector), spawn_positions)

    # Simulation parameters
    time_step = T(1)u"ns"
    max_nsteps = 10000
    convergence_threshold = 1 - 1e-6
    intersect_op = Intersect(mintot = 0)

    n_points = length(inside_detector_idx)
    drift_times = Vector{Int}(undef, n_points)

    @info "Simulating $n_points energy depositions on grid r=0:$grid_step:$radius, z=0:$grid_step:$height at $(angle_deg)°..."

    @threads for i in 1:n_points
        if i % 1000 == 0
            pct = round(100 * i / n_points)
            println("...simulating $i / $n_points ($pct%)")
        end

        pos = find_valid_spawn_position(inside_detector_idx[i], spawn_positions, sim.detector; verbose = false)
        event = SSD.Event([pos], [2039u"keV"])
        simulate!(event, sim, Δt = time_step, max_nsteps = max_nsteps, verbose = false)

        wf = get_electron_and_hole_contribution(event, sim, 1).hole_contribution
        drift_times[i] = extract_drift_time_from_waveform(ustrip(wf.signal), convergence_threshold, intersect_op)
    end

    # Build drift time matrix (r × z layout)
    drift_time_matrix = fill(NaN, length(r_axis), length(z_axis))
    for (i, idx) in enumerate(idx_spawn_positions[inside_detector_idx])
        drift_time_matrix[idx] = drift_times[i]
    end

    # Transpose to z × r layout for output, then extend drift time map and axes
    transposed_map = transpose(drift_time_matrix)  # Now z × r
    r_axis_with_units = collect(r_axis) * u"m"
    z_axis_with_units = collect(z_axis) * u"m"

    # extend_drift_time_map expects (rows × cols) with corresponding (row_axis, col_axis)
    # transposed_map is (z × r), so pass (z_axis, r_axis)
    # Note: col_axis (r) won't be extended into negative values
    extended = extend_drift_time_map(transposed_map, z_axis_with_units, r_axis_with_units)

    ang_str = lpad(string(angle_deg), 3, '0')
    return (;
        :r => extended.col_axis,  # column axis of extended map (r)
        :z => extended.row_axis,  # row axis of extended map (z)
        Symbol("drift_time_$(ang_str)_deg") => extended.drift_map * u"ns"
    )
end
