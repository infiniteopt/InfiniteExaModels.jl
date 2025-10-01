@testset "Test Problem 1" begin
    # Setup model and get ground truth results
    m = InfiniteModel(Ipopt.Optimizer)
    set_silent(m)
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t), x) + 2y(0, 1))
    @constraint(m, ∂(y, t) == sin(y) + z)
    @constraint(m, y + z <= 42, DomainRestrictions(t => [0, 0.5]))
    @constraint(m, y(0, x) == 5)
    optimize!(m)
    obj = objective_value(m)
    yval = value(y)
    zval = value(z)
    # test with ExaTranscriptionBackend
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = 1e-6)
    @test all(isapprox.(yval, value(y), atol = 1e-6))
    @test isapprox.(zval, value(z), atol = 1e-6)
end