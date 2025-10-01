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

@testset "Parameter Function Problem" begin
    # Test custom function for parameter function
    function paramFunc2(t, s, tk)
        if t ≤ 0.5
            return cos(t)*s - tk
        else
            return sin(t)*s + tk
        end
    end

    # Build the InfiniteModel
    ti = 0.2
    model = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
    @infinite_parameter(model, t ∈ [0, 1], num_supports = 5)
    @infinite_parameter(model, s ∈ [2, 3], num_supports = 5)
    @variable(model, 0 ≤ v ≤ 100, Infinite(t))
    @parameter_function(model, pf == sin(t))
    @parameter_function(model, pf2 == (t, s) -> paramFunc2(t, s, ti))
    @constraint(model, c1, v + pf ≤ 100)
    @constraint(model, c2, v*2 + pf*pf2 ≤ 100)
    @constraint(model, c3, v ≥ 0.5*pf2)
    @objective(model, Min, ∫(v*pf, t) + ∫(∫(0.5*v*pf2, t), s))
    set_silent(model)
    optimize!(model)
    tol = 1E-6
    @test isapprox(objective_value(model), 1.907716395171399, atol=tol)
    @test isapprox(
        value(v),
        [
            1.4000000045963952,
            1.3533686294198077,
            1.2163738391341032,
            1.1224581356573253,
            1.362206476152603,
        ],
    )
end