
using Test
using SolidStateDetectors
using LegendSimflow
using LinearAlgebra  # for norm
using PropDicts

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

    meta, xtal, opv = load_detector_metadata(meta, det, opv_val)

    @test isa(xtal, PropDict)
    @test isa(meta, PropDict)

    @test opv_val == opv
end

@testset "setup_hpge_simulation" begin
    meta_path = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "tests", "dummyprod", "inputs"))
    det = "V99000A"
    opv_val = 3000.0
    T = Float64
    refinement_limits = [0.0, 10.0]

    meta, xtal, opv = load_detector_metadata(meta_path, det, opv_val)

    sim = setup_hpge_simulation(meta_path, meta, xtal, opv, T, refinement_limits, threshold = 2000)

    @test isa(sim, Simulation{T})
end
