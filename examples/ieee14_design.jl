using Revise
Revise.includet("../InfiniteExaModels.jl") # allows updates in functions, not data --> eventually since should just be a package
using .InfiniteExaModels

using InfiniteOpt, Distributions, NLPModelsIpopt, Ipopt, Random
Random.seed!(42)

# Set the parameters
θ_nom = [0.; 60.; 10.]
covar = [80. 0 0; 0 80. 0; 0 0 120.]
n_z = 3; n_θ = 3; n_d = 3
c = ones(n_d) / sqrt(n_d)
c_max = 5
U = 10000

# Define the infinite model
im = InfiniteModel(Ipopt.Optimizer)
set_silent(im)
@infinite_parameter(im, θ[i = 1:n_θ] ~ MvNormal(θ_nom, covar), num_supports = 1000)
@variable(im, 0 <= y <= 1, Infinite(θ))
@variable(im, z[1:n_z], Infinite(θ))
@variable(im, d[1:n_d] >= 0)
@objective(im, Max, expect(1 - y, θ))
@constraint(im, f1, -z[1] - 35 - d[1] <= y * U)
@constraint(im, f2, z[1] - 35 - d[1] <= y * U)
@constraint(im, f3, -z[2] - 50 - d[2] <= y * U)
@constraint(im, f4, z[1] - 50 - d[2] <= y * U)
@constraint(im, f5, -z[3] <= y * U)
@constraint(im, f6, z[3] - 100 - d[3] <= y * U)
@constraint(im, h1, z[1] - θ[1] == 0)
@constraint(im, h2, -z[1] -z[2] + z[3] - θ[2] == 0)
@constraint(im, h3, z[2] - θ[3] == 0)
@constraint(im, max_cost, sum(c[i] * d[i] for i = 1:n_d) <= c_max)

# Create the ExaModel and solve both models to compare
@time em, mappings = exa_model(im)
optimize!(im)
result = ipopt(em, print_level = 0)

# Get the answers
ed = [result.solution[mappings.finvar_mappings[v].i] for v in d]
id = value.(d)

# Print a report
println("\n--------------------------------------------")
println("               SUMMARY")
println("--------------------------------------------\n")
println("ExaModel Objective:      ", -result.objective) # change sign for maximization
println("InfiniteModel Objective: ", objective_value(im))
println("\nExaModel d:      ", ed)
println("InfiniteModel d: ", id)
