using Printf

seed = 0
filename = "pglib_opf_case3_lmbd.m"

include("utils.jl")
include("opf.jl")
include("quadrotor.jl")

num_supports_list = [1000,2000,4000,8000,16000]

# for num_supports = num_supports_list
#     quad(filenum_supports = num_supports)
#     opf(filename; seed = seed, num_supports = num_supports)
# end

opf_nvars = []
opf_ncons = []
opf_exa_total_time = []
opf_exa_ad_time = []
opf_jump_total_time = []
opf_jump_ad_time = []

quad_nvars = []
quad_ncons = []
quad_exa_total_time = []
quad_exa_ad_time = []
quad_jump_total_time = []
quad_jump_ad_time = []

for num_supports = num_supports_list
    nvar, ncon, tot, ad = ipopt_stats(joinpath("logs","opf_exa_$(filename)_$(num_supports)_$(seed).log"))
    push!(opf_exa_total_time, tot)
    push!(opf_exa_ad_time, ad)

    nvar, ncon, tot, ad = ipopt_stats(joinpath("logs","opf_jump_$(filename)_$(num_supports)_$(seed).log"))
    push!(opf_jump_total_time, tot)
    push!(opf_jump_ad_time, ad)
    push!(opf_nvars, nvar)
    push!(opf_ncons, ncon)

    nvar, ncon, tot, ad = ipopt_stats(joinpath("logs","quad_exa_$(num_supports).log"))
    push!(quad_exa_total_time, tot)
    push!(quad_exa_ad_time, ad)

    nvar, ncon, tot, ad = ipopt_stats(joinpath("logs","quad_jump_$(num_supports).log"))
    push!(quad_jump_total_time, tot)
    push!(quad_jump_ad_time, ad)

    push!(quad_nvars, nvar)
    push!(quad_ncons, ncon)
end


tbl = join((
"""
$(mod(i,5) == 1 ? "\\hline" : "")
$(varcon(nsup))
& $(varcon(opf_nvars[i]))
& $(varcon(opf_ncons[i]))
& $(fmt(opf_exa_ad_time[i]))
& $(fmt(opf_exa_total_time[i]))
& $(fmt(opf_jump_ad_time[i]))
& $(fmt(opf_jump_total_time[i]))
& $(varcon(nsup))
& $(varcon(quad_nvars[i]))
& $(varcon(quad_ncons[i]))
& $(fmt(quad_exa_ad_time[i]))
& $(fmt(quad_exa_total_time[i]))
& $(fmt(quad_jump_ad_time[i]))
& $(fmt(quad_jump_total_time[i]))
"""
        for (i, nsup) in enumerate(num_supports_list)
), "\\\\\n")


write(
    "results/results.tex",
    replace(
        read("template.tex", String),
        "%% data %%" => replace(
            tbl,
            "_" => "\\_"
        )
    )
)

run(`pdflatex results/results.tex`)
