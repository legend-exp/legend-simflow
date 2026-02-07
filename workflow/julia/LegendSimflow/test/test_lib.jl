
using Test
using SolidStateDetectors
using LegendSimflow
using LinearAlgebra  # for norm

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
    interior = axis[2:end-1]
    @test all(x -> x > 0, interior)
    @test all(x -> x < boundary, interior)

    # Interior spacing is uniform
    steps = diff(interior)
    @test all(isapprox.(steps, steps[1]))

    # Grid step is close to requested value (up to rounding)
    @test isapprox(steps[1], grid_step; rtol=1e-12)
end