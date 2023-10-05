include("InfiniteExaModels.jl")
using InfiniteOpt, Distributions, NLPModelsIpopt, Ipopt, Random
Random.seed!(42)

# Set the data
num_scenarios = 10 # small amount for example
α = [150, 230, 260] # land cost
β = [238, 210, 0]   # purchasing cost
λ = [170, 150, 36]  # selling price
d = [200, 240, 0]   # contract demand
xbar = 500          # total land
wbar3 = 6000        # no upper bound on the other crops
ybar3 = 0           # no upper bound on the other crops
Ξ = [Uniform(0, 5), Uniform(0, 5), Uniform(10, 30)]; # the distributions

# Define the InfiniteModel
im = InfiniteModel(Ipopt.Optimizer)
set_silent(im)
@infinite_parameter(im, ξ[c in 1:3] ~ Ξ[c], num_supports = num_scenarios)
@variables(im, begin
    # 1st stage variables
    0 <= x[1:3] <= xbar
    # 2nd stage variables
    0 <= y[1:3], Infinite(ξ)
    0 <= w[1:3], Infinite(ξ)
end)
@objective(im, Min, α'x + 𝔼(β'y - λ'w, ξ))
@constraints(im, begin
    # capacity constraint
    sum(x) <= xbar
    # balances
    ξ .* x + y - w .>= d
    # crop limits
    w[3] <= wbar3
    y[3] <= ybar3
end)

# Create the ExaModel and solve both models to compare
em, mappings = exa_model(im)
optimize!(im)
result = ipopt(em, print_level = 0)

# Get the answers
ex = [result.solution[mappings.finvar_mappings[v].i] for v in x]
ix = value.(x)

# Print a report
println("\n--------------------------------------------")
println("               SUMMARY")
println("--------------------------------------------\n")
println("ExaModel Objective:      ", result.objective)
println("InfiniteModel Objective: ", objective_value(im))
println("\nExaModel x:      ", ex)
println("InfiniteModel x: ", ix)
