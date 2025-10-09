tol = 1E-6
@testset "Test Problem 1" begin
    # Setup model and get ground truth results
    m = InfiniteModel(Ipopt.Optimizer)
    set_silent(m)
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t), x) + 2y(0, 1))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t, DomainRestrictions(t => [0, 0.5]))
    @constraint(m, ∂(y(0, x), x) == 5)
    optimize!(m)
    obj = objective_value(m)
    yval = value(y)
    zval = value(z)
    # test with ExaTranscriptionBackend
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
    # Test Orthogonal Collocation
    set_derivative_method(t, OrthogonalCollocation(3))
    set_transformation_backend(m, TranscriptionBackend(Ipopt.Optimizer))
    @variable(m, u, Infinite(t))
    constant_over_collocation(u, t)
    set_silent(m)
    optimize!(m)
    obj = objective_value(m)
    yval = value(y)
    zval = value(z)
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
end

@testset "Test Problem 2" begin
    # Setup model and get ground truth results
    m = InfiniteModel(Ipopt.Optimizer)
    set_silent(m)
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t) + 2z, x) + 2y(0, 1))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t)
    @constraint(m, ∂(y(0, x), x) == 5)
    optimize!(m)
    obj = objective_value(m)
    yval = value(y)
    zval = value(z)
    # test with ExaTranscriptionBackend
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
    # Test other objective
    @objective(m, Min, ∫(∫(y^2, t) + 2z^2, x) + 2y(0, 1))
    set_transformation_backend(m, TranscriptionBackend(Ipopt.Optimizer))
    set_silent(m)
    optimize!(m)
    obj = objective_value(m)
    yval = value(y)
    zval = value(z)
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
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

    # Build the InfiniteModel & get ground truth results
    ti = 0.2
    m = InfiniteModel(Ipopt.Optimizer)
    @infinite_parameter(m, t ∈ [0, 1], num_supports = 5)
    @infinite_parameter(m, s ∈ [2, 3], num_supports = 5)
    @variable(m, 0 ≤ v ≤ 100, Infinite(t))
    @variable(m, 0 ≤ z ≤ 100, Infinite(t, s))
    @parameter_function(m, pf == sin(t))
    @parameter_function(m, pf2 == (t, s) -> paramFunc2(t, s, ti))
    @constraint(m, c1, v + pf ≤ 100)
    @constraint(m, c2, v*2 + pf*pf2 ≤ 100)
    @constraint(m, c3, v ≥ 0.5*pf2)
    @constraint(m, c4, z(t, 2.5) + pf2*pf ≤ 40)    # Test semi-infinite variable in constraint
    @constraint(m, c5, v*∫(pf2, s) ≤ 100)   # Test semi-infinite param func in measure
    @objective(m, Min, ∫(v*pf, t) + ∫(∫(0.5*z*pf2, t), s))
    set_silent(m)
    optimize!(m)
    obj = objective_value(m)
    vVal = value(v)
    zVal = value(z)
    # test with ExaTranscriptionBackend
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(vVal, value(v), atol = tol))
    @test all(isapprox.(zVal, value(z), atol = tol))
end