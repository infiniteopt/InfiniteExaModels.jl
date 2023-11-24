function ipopt_stats(fname)
    
    output = read(fname, String)
    tot = parse(Float64, split(split(split(output, "OverallAlgorithm....................:")[2], "wall:")[2], ")")[1])
    ad = parse(Float64,split(split(split(output, "Function Evaluations................:")[2], "wall:")[2], ")")[1])
    nvar = parse(Float64, split(split(output,"Total number of variables............................:")[2], "\n")[1])
    necon = parse(Int, split(split(output,"Total number of equality constraints.................:")[2], "\n")[1])
    nicon = parse(Int, split(split(output,"Total number of inequality constraints...............:")[2], "\n")[1])
    ncon = necon + nicon

    return ncon, ncon, tot, ad
end

function varcon(n)
    if n < 1000
        "$n"
    elseif n < 1000000'
        @sprintf("%5.1fk", n/1000)
    else
        @sprintf("%5.1fm", n/1000000)
    end
end
fmt(t) = @sprintf("%5.2f", t)
efmt(t) = @sprintf("%1.8e", t)
percent(t) = @sprintf("%5.1f", t * 100) * "\\%"
