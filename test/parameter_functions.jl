@testset "Parameter Functions" begin
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