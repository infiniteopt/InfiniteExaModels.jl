using Revise
Revise.includet("./InfiniteExaModels.jl") # allows updates in functions, not data --> eventually since should just be a package
using .InfiniteExaModels

using InfiniteOpt, Distributions, Ipopt

im = InfiniteModel(Ipopt.Optimizer)
@infinite_parameter(im, t in [0, 1], num_supports = 4)
@infinite_parameter(im, ξ[1:2] ~ Uniform(0, 1), num_supports = 3)
@variable(im, y >= 1, Infinite(t), start = 1)
@variable(im, q <= 0.2, Infinite(ξ, t))
@variable(im, z == 5)
@parameter_function(im, b == (t, ξ) -> t * sum(ξ))
@constraint(im, c1, sin(y) + 2z + t == 0)
@constraint(im, c2, q^2 - y + ξ[2] <= 2)
@constraint(im, c3, z^3 + y(0) >= 2)
@constraint(im, c4, 0 <= q(ξ, 0) <= 1) # ExaModel is not detecting this as a linear constraint...
@constraint(im, c5, q(supports(ξ)[:, 2], t)^4 <= 1)
@objective(im, Min, z + ∫(𝔼((b - q)^2, ξ), t))

em, d = exa_model(im);
