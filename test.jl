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
@constraint(im, c3, z^3 + y(0) >= 2)
@constraint(im, c4, 0 <= q(ξ, 0) <= 1) # ExaModel is not detecting this as a linear constraint...
@constraint(im, c5, q(supports(ξ)[:, 2], t)^4 <= 1)
@objective(im, Max, z)

em, d = exa_model(im);
