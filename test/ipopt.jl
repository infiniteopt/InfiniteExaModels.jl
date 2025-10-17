@testset "Ipopt option updates" begin
    # Create base problem with some solver options
    m = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t)
    set_silent(m)
    set_time_limit_sec(m, 120.0)
    @test m.backend.silent == true
    @test m.backend.time_limit == 120.0
    optimize!(m)
    tol = 1e-6
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    println("m.backend.options = $(m.backend.options)")
    @test !isempty(m.backend.prev_options)
    @test m.backend.prev_options == Dict(:print_level => 0, :max_wall_time => 120.0)

    # Update & add new solver options
    set_optimizer_attribute(m, :max_iter, 50)
    set_optimizer_attribute(m, :mu_init, 1e-2)
    set_optimizer_attribute(m, :max_wall_time, 200.0)
    set_optimizer_attribute(m, :tol, 1e-6)  # new option

    # Resolve the same problem
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    @test m.backend.prev_options == Dict(:max_iter => 50, :mu_init => 1e-2, :max_wall_time => 200.0, :print_level => 0, :tol => 1e-6)
end