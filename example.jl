include("InfiniteExaModels.jl")
using InfiniteOpt, Distributions, NLPModelsIpopt, Ipopt

im = InfiniteModel(Ipopt.Optimizer)
@infinite_parameter(im, ξ[1:2] ~ Uniform(0, 1), num_supports = 3)
@variable(im, y >= 0, Infinite(ξ))
@variable(im, z)
@constraint(im, y + z <= 10)
# optimize!(im)

em, d = exa_model(im)
# result = ipopt(em)