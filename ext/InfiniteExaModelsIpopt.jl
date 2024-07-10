module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt 

function InfiniteExaModels.translate_option(p::Pair, ::Type{NLPModelsIpopt.IpoptSolver})
    if p[1] == :verbose
        return :print_level => isone(p[2]) ? 5 : 0
    else
        return p
    end
end

end