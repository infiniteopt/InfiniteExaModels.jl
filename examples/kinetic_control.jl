using InfiniteExaModels

using InfiniteOpt, NLPModelsIpopt, Ipopt

# Define parameters
A = [3.6362e6, 2.5212e16, 190.6879, 8.7409e24]
Ea = [10000, 25000, 5000, 40000]
R = 1.987
T_lower = 273. + 40
T_upper = 273. + 60
c0 = [1., 0., 0.]
Tr = 273 .+ [30, 40, 50, 70]
kr = [A[j] * exp(-Ea[j] / R / Tr[j]) for j in eachindex(A)]
tf = 3

# Define the model
im = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
@infinite_parameter(im, t ∈ [0, tf], num_supports = 100, derivative_method = OrthogonalCollocation(4))
add_supports(t, [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1]) # add initial supports with finer resolution
@variable(im, 0 ≤ c[i = 1:3] ≤ 1, Infinite(t), start = c0[i])
@variable(im, T_lower ≤ T ≤ T_upper, Infinite(t), start = T_upper)
@objective(im, Max, c[2](tf))
@constraint(im, [i = 1:3], c[i](0) == c0[i])
@expression(im, k[j = 1:4], kr[j] * exp(Ea[j] / R * (1 / Tr[j] -  1 / T))) # compute k relative to kr (k at a reference temp) for better scaling
@expression(im, r1, c[1] * k[1] - c[2] * k[2])
@expression(im, r2, c[1] * k[3] - c[3] * k[4])
@constraint(im, b1, ∂(c[1], t) == -r1 - r2)
@constraint(im, ∂(c[2], t) == r1)
@constraint(im, ∂(c[3], t) == r2)
constant_over_collocation(T, t)

# Solve
set_attribute(im, "print_timing_statistics", "yes")
optimize!(im)
