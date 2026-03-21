
using Test
using SolidStateDetectors
using LegendSimflow
using LinearAlgebra  # for norm
using PropDicts
using RadiationDetectorDSP
using Unitful
using LegendDataManagement

@testset "build_simulation_grid_axis" begin
    T = Float64
    boundary = 10.0
    grid_step = 1.0

    axis = build_simulation_grid_axis(boundary, grid_step)

    tol = SolidStateDetectors.ConstructiveSolidGeometry.csg_default_tol(T)
    offset = 2 * tol

    @test isa(axis, Vector{T})

    # Check extended boundary points
    @test axis[1] ≈ -offset
    @test axis[end] ≈ boundary + offset

    # Interior points lie strictly within (0, boundary)
    interior = axis[2:(end - 1)]
    @test all(x -> x > 0, interior)
    @test all(x -> x < boundary, interior)

    # Interior spacing is uniform
    steps = diff(interior)
    @test all(isapprox.(steps, steps[1]))

    # Grid step is close to requested value (up to rounding)
    @test isapprox(steps[1], grid_step; rtol = 1e-12)
end

@testset "load_detector_metadata" begin

    meta = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "tests", "dummyprod", "inputs"))

    det = "V99000A"
    opv_val = 3000.0

    meta_dict, xtal, opv = load_detector_metadata(meta, det, opv_val)

    @test isa(xtal, PropDict)
    @test isa(meta_dict, PropDict)

    @test opv_val == opv

    # test with no opv
    meta_dict, xtal, opv = load_detector_metadata(meta, det)

    @test isa(xtal, PropDict)
    @test isa(meta_dict, PropDict)

    # opv should be read from metadata and match provided value
    @test opv == Float32(meta_dict.characterization.l200_site.recommended_voltage_in_V)

end

@testset "map_generation" begin
    meta_path = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "tests", "dummyprod", "inputs"))
    det = "V99000A"
    opv_val = 3000.0
    T = Float64
    refinement_limits = [0.2, 0.1, 0.05, 0.02]

    meta, xtal, opv = load_detector_metadata(meta_path, det, opv_val)

    sim = setup_hpge_simulation(meta_path, meta, xtal, opv, T, refinement_limits, threshold = 20000)

    @test isa(sim, Simulation{T})

    # test valid spawn pos
    # Candidate positions include one valid and one invalid
    spawn_positions =
        [CartesianPoint(T(20/1000.0), T(0.0), T(-5/1000.0)), CartesianPoint(T(20/1000.0), T(0.5/1000.0), T(5/1000.0))]

    # Test that the valid position is returned
    pos = find_valid_spawn_position(2, spawn_positions, sim.detector; verbose = false)
    @test pos == spawn_positions[2]

    pos = find_valid_spawn_position(1, spawn_positions, sim.detector; verbose = false)
    @test pos == spawn_positions[2]

    dt_map = compute_drift_time_map(sim, meta, T, 0.0, 10/1000.0, 0)

    @test hasproperty(dt_map, :drift_time_000_deg)
    @test hasproperty(dt_map, :r)
    @test hasproperty(dt_map, :z)

    @test size(dt_map.drift_time_000_deg) == (length(dt_map.z), length(dt_map.r))

    wf_map = compute_ideal_pulse_shape_lib(sim, meta, T, 0.0, false, 10/1000.0)

    @test hasproperty(wf_map, :waveform_000_deg)
    @test hasproperty(wf_map, :r)
    @test hasproperty(wf_map, :z)
    @test hasproperty(wf_map, :dt)
    @test wf_map.dt == 1.0 .* u"ns"

    @test size(wf_map.waveform_000_deg) == (5000, length(wf_map.z), length(wf_map.r))

end

@testset "extract_drift_time_from_waveform" begin
    convergence_threshold = 1 - 1e-6
    intersect_op = Intersect(mintot = 0)

    # Waveform that ramps from 0 to 1 then holds at plateau
    n = 100
    wf_positive = vcat(collect(1:50) ./ 50.0, ones(50))

    dt = extract_drift_time_from_waveform(wf_positive, convergence_threshold, intersect_op)

    @test isa(dt, Real)
    @test 1 <= dt <= n

    # Negative waveform: negation is handled internally and should give the same result
    wf_negative = -wf_positive
    dt_neg = extract_drift_time_from_waveform(wf_negative, convergence_threshold, intersect_op)

    @test isa(dt_neg, Real)
    @test dt_neg == dt
end

@testset "extend_drift_time_map" begin
    # 3×3 drift map with NaN at corners; interior values are finite
    drift_map = [
        NaN 2.0 NaN
        1.0 3.0 4.0
        NaN 5.0 NaN
    ]
    row_axis = [0.0, 1.0, 2.0] * u"m"
    col_axis = [0.0, 1.0, 2.0] * u"m"

    result = extend_drift_time_map(drift_map, row_axis, col_axis; layers = 1)

    # Return value is a NamedTuple with the expected fields
    @test result isa NamedTuple
    @test hasproperty(result, :drift_map)
    @test hasproperty(result, :row_axis)
    @test hasproperty(result, :col_axis)

    # With layers=1: rows extended on both sides (+2 total), cols on high side only (+1)
    @test size(result.drift_map) == (5, 4)
    @test length(result.row_axis) == 5
    @test length(result.col_axis) == 4

    # Original non-NaN values are preserved (original data placed at rows 2:4, cols 1:3)
    @test result.drift_map[3, 2] ≈ 3.0  # original centre value

    # NaN corner at original [1,1] (work_map[2,1]) has non-NaN neighbours and is filled in
    @test !isnan(result.drift_map[2, 1])
end
