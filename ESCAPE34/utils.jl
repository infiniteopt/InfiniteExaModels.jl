using Printf

function ipopt_stats(fname)
    
    output = read(fname, String)
    tot = parse(Float64, split(split(split(output, "OverallAlgorithm....................:")[2], "wall:")[2], ")")[1])
    ad = parse(Float64,split(split(split(output, "Function Evaluations................:")[2], "wall:")[2], ")")[1])
    nvar = parse(Int, split(split(output,"Total number of variables............................:")[2], "\n")[1])
    necon = parse(Int, split(split(output,"Total number of equality constraints.................:")[2], "\n")[1])
    nicon = parse(Int, split(split(output,"Total number of inequality constraints...............:")[2], "\n")[1])
    ncon = necon + nicon

    return nvar, ncon, tot, ad
end

function madnlp_stats(fname)
    
    output = read(fname, String)
    tot = -1.0
    ad = -1.0
    try
        tot = parse(Float64, split(split(output,"Total wall-clock secs                                       =")[2], "\n")[1])
        ad = parse(Float64, split(split(output,"Total wall-clock secs in NLP function evaluations           =")[2], "\n")[1])
    catch
        
    end
    nvar = parse(Int, split(split(output,"Total number of variables............................:")[2], "\n")[1])
    necon = parse(Int, split(split(output,"Total number of equality constraints.................:")[2], "\n")[1])
    nicon = parse(Int, split(split(output,"Total number of inequality constraints...............:")[2], "\n")[1])
    ncon = necon + nicon

    return nvar, ncon, tot, ad
end
