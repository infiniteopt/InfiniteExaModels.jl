include("opf.jl")
include("quadrotor.jl")

for num_supports = [100, 1000, 10000, 100000]
    quad(num_supports = num_supports)
    opf(num_supports = num_supports)
end
