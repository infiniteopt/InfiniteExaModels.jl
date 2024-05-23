using InfiniteExaModels

using InfiniteOpt, NLPModelsIpopt, Ipopt

# Data
xw = [1 4 6 1; 1 3 0 1] # positions
tw = [0, 25, 50, 60];    # times
dmethod = OrthogonalCollocation(4)
dmethod = FiniteDifference(Backward())

# Define the InfiniteModel
im = InfiniteModel(Ipopt.Optimizer)
set_silent(im)
@infinite_parameter(im, t in [0, 60], num_supports = 101, derivative_method = dmethod)
@variables(im, begin
    # state variables
    x[1:2], Infinite(t)
    v[1:2], Infinite(t)
    # control variables
    u[1:2], Infinite(t), (start = 0)
end)
@objective(im, Min, ∫(u[1]^2 + u[2]^2, t))
@constraint(im, [i = 1:2], v[i](0) == 0)
@constraint(im, [i = 1:2], ∂(x[i], t) == v[i])
@constraint(im, [i = 1:2], ∂(v[i], t) == u[i])
@constraint(im, [i = 1:2, j = eachindex(tw)], x[i](tw[j]) == xw[i, j])
constant_over_collocation.(u, t) # needed for collocation

# Create the ExaModel and solve both models to compare
@time em, mappings = exa_model(im)
optimize!(im)
result = ipopt(em, print_level = 0)

# Print a report
println("\n--------------------------------------------")
println("               SUMMARY")
println("--------------------------------------------\n")
println("ExaModel Objective:      ", result.objective)
println("InfiniteModel Objective: ", objective_value(im))
