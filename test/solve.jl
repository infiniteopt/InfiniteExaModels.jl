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
    dyval = value(∂(y, t))
    # test with ExaTranscriptionBackend
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
    @test isapprox(dyval, value(∂(y, t)), atol=tol)
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
    dyval = value(∂(y, t))
    @test set_transformation_backend(m, ExaTranscriptionBackend(IpoptSolver)) isa Nothing
    @test set_silent(m) isa Nothing
    @test optimize!(m).status == :first_order 
    @test isapprox(obj, objective_value(m), atol = tol)
    @test all(isapprox.(yval, value(y), atol = tol))
    @test isapprox.(zval, value(z), atol = tol)
    @test isapprox(dyval, value(∂(y, t)), atol=tol)
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
    # Test other objectives
    int = ∫(y^2, t)
    objs = [
        ∫(int + 2z^2, x) + 2y(0, 1),
        ∫(int + sin(z^2), x),
        ∫(int * cos(z), x),
        ∫(z * (int + z^3), x)
    ]
    for obj in objs
        @objective(m, Min, obj)
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
    @constraint(m, c3, v ≥ 0.2*pf2)
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
    @test all(isapprox.(vVal, value(v), atol = 1e-4))
    @test all(isapprox.(zVal, value(z), atol = tol))
end

@testset "Parameter updates" begin
    m = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
    @infinite_parameter(m, t in [0, 1], num_supports = 3)
    @finite_parameter(m, p1 == 100.0)
    @finite_parameter(m, p2 == 1.0)
    @variable(m, x[1:2], Infinite(t))
    @objective(m, Min, p1 * ∫((x[2] - x[1]^2)^2, t) + ∫((p2 - x[1])^2, t))
    @constraint(m, [i in 1:2], x[i] ≤ [0.5, 3.0][i])
    @constraint(m, x[1] * x[2] ≥ 1.0)
    @constraint(m, x[1] + x[2]^2 ≥ 0.0)
    set_silent(m)
    optimize!(m)
    @test isapprox(objective_value(m), 306.4999755050365, atol=tol)
    @test value(p1) == 100.0
    @test value(p2) == 1.0
    # Update the parameters
    set_parameter_value(p1, 90.0)
    set_parameter_value(p2, 1.3)
    # Solve again
    optimize!(m)
    @test isapprox(objective_value(m), 276.26497794903645, atol=tol)
    @test value(p1) == 90.0
    @test value(p2) == 1.3
    # Add a new finite parameter & try to update
    @finite_parameter(m, p3 == 43.0)
    @constraint(m, x[1]^2 + x[2]^2 ≤ p3)
    set_parameter_value(p3, 50.0)
    @test transformation_backend_ready(m) == false
end

@testset "Parameter function updates" begin
    # Define function for updating parameter function later
    function oldpf2(t, s)
        return sin(t)*s + 0.2
    end
    function newpf2(t, s)
        return sin(t)*s + 0.8
    end
    # Build the InfiniteModel
    m = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
    @infinite_parameter(m, t ∈ [0, 1], num_supports = 3)
    @infinite_parameter(m, s ∈ [2, 3], num_supports = 3)
    @variable(m, 0 ≤ v ≤ 100, Infinite(t))
    @variable(m, 0 ≤ z ≤ 100, Infinite(t, s))
    @parameter_function(m, pf1 == sin(t))
    @parameter_function(m, pf2 == oldpf2(t, s))
    @constraint(m, c1, v + pf1 ≤ 100)
    @constraint(m, c2, v*2 + pf1*pf2 ≤ 100)
    @constraint(m, c3, v ≥ 0.5*pf2)
    @constraint(m, c4, z(t, 2.5) + pf2*pf1 ≤ 40)
    @objective(m, Min, ∫(v*pf1, t) + ∫(∫(0.5*z*pf2, t), s))
    set_silent(m)
    optimize!(m)
    @test isapprox(objective_value(m), 0.48292223509341475, atol=tol)
    param1 = transformation_variable(pf1)
    param2 = transformation_variable(pf2)
    expectedpf1 = sin.([0.0, 0.5, 1.0])
    expectedpf2 = [0.2, 1.158851077208406, 1.882941969615793, 0.2, 1.3985638465105075, 2.3036774620197416, 0.2, 1.638276615812609, 2.7244129544236895]
    em = InfiniteOpt.transformation_model(m)
    @test em.θ[param1.offset+1:param1.offset+param1.length] == expectedpf1
    @test em.θ[param2.offset+1:param2.offset+param2.length] == expectedpf2
    # Test value queries
    @test value(pf1) == expectedpf1
    @test reshape(value(pf2), 9) == expectedpf2
    # Update parameter functions
    set_parameter_value(pf1, cos)
    set_parameter_value(pf2, newpf2)
    expectedpf1 = cos.([0.0, 0.5, 1.0])
    expectedpf2 = [0.8, 1.758851077208406, 2.4829419696157933, 0.8, 1.9985638465105076, 2.9036774620197416, 0.8, 2.238276615812609, 3.324412954423689]
    @test em.θ[param1.offset+1:param1.offset+param1.length] == expectedpf1
    @test em.θ[param2.offset+1:param2.offset+param2.length] == expectedpf2
    optimize!(m)
    @test isapprox(objective_value(m), 0.8155916466182952, atol=tol)
    @test value(pf1) == expectedpf1
    @test reshape(value(pf2), 9) == expectedpf2
end