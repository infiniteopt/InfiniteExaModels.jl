# Store the mappings between JuMP's operator names and the functions supported by ExaModels
const _op_mappings = Dict{Symbol, Function}(
    :+ => +,
    :- => -,
    :* => *,
    :/ => /,
    :^ => ^,
    :inv => inv,
    :sqrt => sqrt,
    :cbrt => cbrt,
    :abs => abs,
    :abs2 => abs2,
    :exp => exp,
    :exp2 => exp2,
    :log => log,
    :log2 => log2,
    :log10 => log10,
    :log1p => log1p,
    :sin => sin,
    :cos => cos,
    :tan => tan,
    :asin => asin,
    :acos => acos,
    :csc => csc,
    :sec => sec,
    :cot => cot,
    :atan => atan,
    :acot => acot,
    :sind => sind,
    :cosd => cosd,
    :tand => tand,
    :cscd => cscd,
    :secd => secd,
    :cotd => cotd,
    :atand => atand,
    :acotd => acotd,
    :sinh => sinh,
    :cosh => cosh,
    :tanh => tanh,
    :csch => csc,
    :sech => sech,
    :coth => coth,
    :atanh => atanh,
    :acoth => acoth,
    # TODO add the remaining JuMP operators
)

# Map a nonlinear function symbol to the underlying function
function _nl_op(s::Symbol)
    if !haskey(_op_mappings, s)
        error("`InfiniteExaModel`s does not support the nonlinear operator `$s`.")
    end
    return _op_mappings[s]
end