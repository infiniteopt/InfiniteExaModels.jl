tol = 1e-6
@testset "MadNLP option updates 1" begin
    # Create base problem with some solver options + silenced
    m = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver))
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t)
    set_silent(m)
    etb = m.backend
    set_time_limit_sec(m, 120.0)
    @test etb.silent == true
    @test etb.time_limit == 120.0
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    @test etb.solver.opt.max_iter == 3000 # MadNLP default
    @test etb.solver.opt.mu_init == 1e-1  # MadNLP default
    @test etb.solver.opt.max_wall_time == etb.time_limit
    @test etb.solver.logger.print_level == MadNLP.ERROR
    @test !isempty(etb.prev_options)
    @test etb.options == Dict(:solver => MadNLPSolver)
    @test etb.prev_options == Dict(:print_level => MadNLP.ERROR, :max_wall_time => 120.0)
    @test !isnothing(etb.results)

    # Update & add new solver options
    unset_silent(m) # Turn off silent mode
    set_time_limit_sec(m, 200.0)   # Change time time_limit
    set_optimizer_attribute(m, :max_iter, 50)
    set_optimizer_attribute(m, :mu_init, 1e-2)
    set_optimizer_attribute(m, :tol, 1e-6)  # new option
    @test etb.silent == false

    # Ensure that the results from before weren't wiped after changing options
    @test !isnothing(etb.results)

    # Resolve the same problem
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    @test etb.solver.opt.max_iter == 50
    @test etb.solver.opt.mu_init == 1e-2
    @test etb.options ==
        Dict(
            :solver   => MadNLPSolver,
            :max_iter => 50,
            :mu_init  => 1e-2,
            :tol      => 1e-6,
        )
    @test etb.prev_options == 
        Dict(
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :print_level => MadNLP.INFO,
            :max_wall_time => 200.0,
        )
    @test etb.solver.logger.print_level == MadNLP.INFO
end

@testset "MadNLP option updates 2" begin
    # Create base problem with some solver options + unsilenced
    m = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver))
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t)
    etb = m.backend
    set_time_limit_sec(m, 120.0)
    @test etb.time_limit == 120.0
    set_optimizer_attribute(m, :max_iter, 50)
    set_optimizer_attribute(m, :mu_init, 1e-2)
    set_optimizer_attribute(m, :tol, 1e-6)
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    @test etb.solver.opt.max_iter == 50 # MadNLP default
    @test etb.solver.opt.mu_init == 1e-2  # MadNLP default
    @test etb.solver.opt.max_wall_time == etb.time_limit
    @test !isempty(etb.prev_options)
    @test etb.options == 
        Dict(
            :solver => MadNLPSolver,
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
        )
    @test etb.prev_options == 
        Dict(
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :max_wall_time => 120.0
        )
    @test !isnothing(etb.results)

    # Change the print level & unset time limit
    set_optimizer_attribute(m, :print_level, MadNLP.WARN)
    unset_time_limit_sec(m)
    @test isnan(etb.time_limit)

    # Solve again
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    prev = etb.prev_options
    @test length(keys(prev)) == 1
    @test prev[:print_level] == MadNLP.WARN
    @test etb.solver.logger.print_level == MadNLP.WARN
    @test etb.options ==
        Dict(
            :solver => MadNLPSolver,
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :print_level => MadNLP.WARN,
        )

    # Set silent mode and time limit
    set_silent(m)
    set_time_limit_sec(m, 150.0)
    @test etb.silent == true
    @test etb.time_limit == 150.0
    optimize!(m)
    @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
    @test etb.solver.logger.print_level == MadNLP.ERROR
    @test etb.options ==
        Dict(
            :solver => MadNLPSolver,
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :print_level => MadNLP.WARN,
        )
    @test etb.prev_options ==
        Dict(
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :print_level => MadNLP.ERROR,
            :max_wall_time => 150.0
        )
    # Restore previous print level after unsetting silent
    unset_silent(m)
    optimize!(m)
    @test etb.options ==
        Dict(
            :solver => MadNLPSolver,
            :max_iter => 50,
            :mu_init => 1e-2,
            :tol => 1e-6,
            :print_level => MadNLP.WARN,
        )
    @test etb.prev_options ==
        Dict(
            :print_level => MadNLP.WARN,
        )
end