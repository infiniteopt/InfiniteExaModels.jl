@testset "MadNLP option updates" begin
    # Create base problem with some solver options
    m = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver))
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
    @test m.backend.solver.opt.max_iter == 3000 # MadNLP default
    @test m.backend.solver.opt.mu_init == 1e-1  # MadNLP default
    @test m.backend.solver.opt.max_wall_time == m.backend.time_limit
    @test !isempty(m.backend.prev_options)
    @test m.backend.prev_options == Dict(:print_level => MadNLP.ERROR, :max_wall_time => 120.0)

    # Update & add new solver options
    unset_silent(m) # Turn off silent mode
    set_time_limit_sec(m, 200.0)   # Change time time_limit
    set_optimizer_attribute(m, :max_iter, 50)
    set_optimizer_attribute(m, :mu_init, 1e-2)
    set_optimizer_attribute(m, :tol, 1e-6)  # new option
    @test m.backend.silent == false

    # Resolve the same problem
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    @test m.backend.solver.opt.max_iter == 50
    @test m.backend.solver.opt.mu_init == 1e-2
    @test m.backend.prev_options == Dict(:max_iter => 50, :mu_init => 1e-2, :max_wall_time => 200.0, :print_level => MadNLP.INFO, :tol => 1e-6)

    # Change the print level & unset time limit
    set_optimizer_attribute(m, :print_level, MadNLP.WARN)
    unset_time_limit_sec(m)
    @test isnan(m.backend.time_limit)

    # Solve one more time
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    prev = m.backend.prev_options
    @test length(keys(prev)) == 5
    @test prev[:max_iter] == 50
    @test prev[:mu_init] == 1e-2
    @test isnan(prev[:max_wall_time])
    @test prev[:print_level] == MadNLP.WARN
    @test prev[:tol] == 1e-6
end