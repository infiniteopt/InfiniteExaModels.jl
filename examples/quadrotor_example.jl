using InfiniteExaModels
using InfiniteOpt, NLPModelsIpopt, Ipopt

function main()

    # Data
    n = 9
    p = 4
    nd = 9

    T = 60
    N = 101

    dmethod = OrthogonalCollocation(3)
    # dmethod = FiniteDifference(Backward())
    
    # Define the InfiniteModel
    im = InfiniteModel(Ipopt.Optimizer)

    @infinite_parameter(im, t in [0, T], num_supports = N, derivative_method = dmethod)
    
    @parameter_function(im, d1 == t -> sin(2 * pi * t/T))
    @parameter_function(im, d3 == t -> 2 * sin(4 * pi * t/T))
    @parameter_function(im, d5 == t -> 2 * (t/T))
    
    @variables(
        im,
        begin
            # state variables
            x[1:n], Infinite(t)
            # control variables
            u[1:p], Infinite(t), (start = 0)
        end
    )
    constant_over_collocation.(u, t)
    @objective(
        im, Min, ∫(
            (x[1] - d1)^2 + (x[3] - d3)^2 + (x[5] - d5)^2 + x[7]^2 + x[8]^2 + x[9]^2
            + 0.1 * (u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2),
            t
        )
    )
    @constraint(im, [i = 1:n], x[i](0) == 0)
    @constraint(
        im, ∂(x[1], t) ==
            x[2]
    )
    @constraint(
        im, ∂(x[2], t) ==
            u[1] * cos(x[7]) * sin(x[8]) * cos(x[9]) + u[1] * sin(x[7]) * sin(x[9])
    )
    @constraint(
        im, ∂(x[3], t) ==
            x[4]            
    )
    @constraint(
        im, ∂(x[4], t) ==
            u[1] * cos(x[7]) * sin(x[8]) * sin(x[9]) - u[1] * sin(x[7]) * cos(x[9])            
    )
    @constraint(
        im, ∂(x[5], t) ==
            x[6]            
    )
    @constraint(
        im, ∂(x[6], t) ==
            u[1] * cos(x[7]) * cos(x[8]) - 9.8            
    )
    @constraint(
        im, ∂(x[7], t) ==
            u[2] * cos(x[7]) / cos(x[8]) + u[3] * sin(x[7]) / cos(x[8])            
    )
    @constraint(
        im, ∂(x[8], t) ==
            -u[2] * sin(x[7]) + u[3] * cos(x[7])
    )
    @constraint(
        im, ∂(x[9], t) ==
            u[2] * cos(x[7]) * tan(x[8]) + u[3] * sin(x[7]) * tan(x[8]) + u[4]
    )

    # Create the ExaModel and solve both models to compare
    @time em, mappings = exa_model(im)
    set_optimizer_attribute(im, "linear_solver", "ma27")
    set_optimizer_attribute(im, "print_timing_statistics", "yes")
    optimize!(im)
    ipopt(em; print_level = 0, linear_solver="ma27")
    result = ipopt(
        em;
        linear_solver="ma27",
        print_timing_statistics = "yes",
    )

    # Print a report
    println("\n--------------------------------------------")
    println("               SUMMARY")
    println("--------------------------------------------\n")
    println("ExaModel Objective:      ", result.objective)
    println("InfiniteModel Objective: ", objective_value(im))

end


main()
