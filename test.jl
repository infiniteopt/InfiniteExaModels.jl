include("InfiniteExaModels.jl")
using InfiniteOpt, Distributions

im = InfiniteModel()
@infinite_parameter(im, t in [0, 1], num_supports = 4)
@infinite_parameter(im, ξ[1:2] ~ Uniform(0, 1), num_supports = 3)
@variable(im, y >= 1, Infinite(t), start = 1)
@variable(im, q <= 2, Infinite(ξ, t))
@variable(im, z == 5)
@constraint(im, c1, sin(y) + 2z + t == 0)
@constraint(im, c2, q^2 - y + ξ[2] <= 2)
@constraint(im, c3, z^3 >= 2)

em, d = exa_model(im)
