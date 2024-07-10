using InfiniteExaModels
using InfiniteOpt

function quad(; num_supports = 100, backend = nothing)
    # Data
    n = 9
    p = 4
    T = 60
    
    # Define the InfiniteModel
    im = InfiniteModel(backend)

    @infinite_parameter(im, t in [0, T], num_supports = num_supports, 
                        derivative_method = OrthogonalCollocation(3))
    
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
    constant_over_collocation.(u, t)

    return im
end
