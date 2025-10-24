tol = 1e-6
# @testset "Ipopt option updates 1" begin
#     # Create base problem with some solver options + silenced
#     m = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
#     @infinite_parameter(m, t in [0, 1], num_supports = 5)
#     @infinite_parameter(m, x in [-1, 1], num_supports = 5)
#     @variable(m, y >= 0, Infinite(t, x))
#     @variable(m, z, start = 10)
#     @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
#     @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
#     @constraint(m, y + z <= 42 + t)
#     set_silent(m)
#     etb = m.backend
#     set_time_limit_sec(m, 120.0)
#     @test etb.silent == true
#     @test etb.time_limit == 120.0
#     @test (output = @capture_out results = optimize!(m)) == ""  
#     @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
#     @test !isempty(etb.prev_options)
#     @test etb.options == Dict(:solver => IpoptSolver)
#     @test etb.prev_options == Dict(:print_level => 0, :max_wall_time => 120.0)
#     @test !isnothing(etb.results)

#     # Update & add new solver options
#     unset_silent(m) # Turn off silent mode
#     set_time_limit_sec(m, 200.0)   # Change time time_limit
#     set_optimizer_attribute(m, :max_iter, 50)
#     set_optimizer_attribute(m, :mu_init, 1e-2)
#     set_optimizer_attribute(m, :tol, 1e-6)  # new option
#     @test etb.silent == false

#     # Check that previous results weren't wiped after changing options
#     @test !isnothing(etb.results)

#     # Resolve the same problem
#     output = @capture_out results = optimize!(m)
#     @test occursin("This is Ipopt version", output)  
#     @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
#     @test etb.options ==
#         Dict(
#             :solver   => IpoptSolver,
#             :max_iter => 50,
#             :mu_init  => 1e-2,
#             :tol      => 1e-6,
#         )
#     @test etb.prev_options ==
#         Dict(
#             :max_iter      => 50,
#             :mu_init       => 1e-2,
#             :tol           => 1e-6,
#             :print_level   => 5,
#             :max_wall_time => 200.0,
#         )
# end

# @testset "Ipopt option updates 2" begin
#     # Create base problem with some solver options + unsilenced
#     m = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
#     @infinite_parameter(m, t in [0, 1], num_supports = 5)
#     @infinite_parameter(m, x in [-1, 1], num_supports = 5)
#     @variable(m, y >= 0, Infinite(t, x))
#     @variable(m, z, start = 10)
#     @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
#     @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
#     @constraint(m, y + z <= 42 + t)
#     etb = m.backend
#     set_time_limit_sec(m, 120.0)
#     @test etb.time_limit == 120.0
#     set_optimizer_attribute(m, :max_iter, 50)
#     set_optimizer_attribute(m, :mu_init, 1e-2)
#     set_optimizer_attribute(m, :tol, 1e-6)
#     output = @capture_out results = optimize!(m)
#     @test occursin("This is Ipopt version", output)
#     @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
#     @test !isempty(etb.prev_options)
#     @test etb.options == 
#         Dict(
#             :solver => IpoptSolver,
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#         )
#     @test etb.prev_options == 
#         Dict(
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#             :max_wall_time => 120.0
#         )
#     @test !isnothing(etb.results)

#     # Change the print level & unset time limit
#     set_optimizer_attribute(m, :print_level, 3)
#     unset_time_limit_sec(m)
#     @test isnan(etb.time_limit)

#     # Solve again
#     output = @capture_out results = optimize!(m)
#     @test occursin("Total number of variables", output)
#     @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
#     prev = etb.prev_options
#     @test prev == Dict(:print_level => 3, :max_wall_time => 1.0e20)
#     @test etb.options ==
#         Dict(
#             :solver => IpoptSolver,
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#             :print_level => 3,
#         )

#     # Set silent mode and time limit
#     set_silent(m)
#     set_time_limit_sec(m, 150.0)
#     @test etb.silent == true
#     @test etb.time_limit == 150.0
#     @test (output = @capture_out optimize!(m)) == ""
#     @test isapprox(objective_value(m), -1.2784599867885884e+01, atol=tol)
#     @test etb.options ==
#         Dict(
#             :solver => IpoptSolver,
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#             :print_level => 3,
#         )
#     @test etb.prev_options ==
#         Dict(
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#             :max_wall_time => 150.0,
#             :print_level => 0,
#         )

#     # Restore previous print level after unsetting silent
#     unset_silent(m)
#     output = @capture_out result = optimize!(m)
#     @test occursin("Total number of variables", output)
#     @test etb.options ==
#         Dict(
#             :solver => IpoptSolver,
#             :max_iter => 50,
#             :mu_init => 1e-2,
#             :tol => 1e-6,
#             :print_level => 3,
#         )
#     @test etb.prev_options ==
#         Dict(
#             :print_level => 3,
#         )
# end

@testset "Ipopt warmstarts" begin
    # Solve base problem
    m = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 5)
    @variable(m, y >= 0, Infinite(t, x))
    @variable(m, z, start = 10)
    @objective(m, Min, ∫(∫(y^2, t) + 2z, x))
    @constraint(m, ∂(y, t) == sin(y) + z + 1.2)
    @constraint(m, y + z <= 42 + t)
    output = @capture_out result = optimize!(m)
    result1 = m.backend.results
    @test occursin("This is Ipopt version", output)
    @test occursin("Number of Iterations....: 8", output)
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    expected = zeros(51)
    expected[1] = 10.0
    @test NLPModels.get_x0(m.backend.model) == expected
    warmstart_backend_start_values(m)
    @test m.backend.options[:x0] == result1.solution
    @test m.backend.options[:y0] == result1.multipliers
    @test m.backend.options[:zL0] == result1.multipliers_L
    @test m.backend.options[:zU0] == result1.multipliers_U
    set_optimizer_attribute(m, "print_user_options", "yes")
    output = @capture_out result = optimize!(m)
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    # Should converge in 4 iterations if warmstarted
    @test occursin("Number of Iterations....: 4", output)
    @test occursin("warm_start_init_point = yes", output)
    # Reset the solution & turn off warmstarting
    m.backend.results = nothing
    set_optimizer_attribute(m, "warm_start_init_point", "no")
    @test_logs (
        :warn,
        "No previous solution values found to warmstart the backend."
        ) warmstart_backend_start_values(m)
    output = @capture_out result = optimize!(m)
    @test m.backend.options[:warm_start_init_point] == "no"
    @test !haskey(m.backend.options, :x0)
    @test !haskey(m.backend.options, :y0)
    @test !haskey(m.backend.options, :zL0)
    @test !haskey(m.backend.options, :zU0)
    @test isapprox(objective_value(m), -1.2784599900757165e+01, atol=tol)
    @test occursin("Number of Iterations....: 8", output)
end