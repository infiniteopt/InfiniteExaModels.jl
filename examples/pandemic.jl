using InfiniteExaModels

using InfiniteOpt, Distributions, NLPModelsIpopt, Ipopt, Random
Random.seed!(42)

# Set the parameters
γ = 0.303
β = 0.727
N = 1e5
extra_ts = [0.001, 0.002, 0.004, 0.008, 0.02, 0.04, 0.08, 0.2, 0.4, 0.8]
dmethod = OrthogonalCollocation(3)
# dmethod = FiniteDifference(Backward())

# Create the infinite model
im = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
set_silent(im)
@infinite_parameter(im, t ∈ [0, 200], num_supports = 101, derivative_method = dmethod)
@infinite_parameter(im, ξ ~ Uniform(0.1, 0.6), num_supports = 4)
add_supports(t, extra_ts)
@variable(im, s ≥ 0, Infinite(t, ξ))
@variable(im, e ≥ 0, Infinite(t, ξ))
@variable(im, i ≥ 0, Infinite(t, ξ))
@variable(im, r ≥ 0, Infinite(t, ξ))
@variable(im, 0 ≤ u ≤ 0.8, Infinite(t), start = 0.2)
@objective(im, Min, ∫(u, t))
@constraint(im, s(0, ξ) == 1 - 1 / N)
@constraint(im, e(0, ξ) == 1 / N)
@constraint(im, i(0, ξ) == 0)
@constraint(im, r(0, ξ) == 0)
@constraint(im, s_constr, ∂(s, t) == -(1 - u) * β * s * i)
@constraint(im, e_constr, ∂(e, t) == (1 - u) * β * s * i - ξ * e)
@constraint(im, i_constr, ∂(i, t) == ξ * e - γ * i)
@constraint(im, r_constr, ∂(r, t) == γ * i)
@constraint(im, imax_constr, i ≤ 0.02)
constant_over_collocation(u, t) # needed for collocation

# Solve
optimize!(im)
