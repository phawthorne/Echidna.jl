module Echidna

include("var_types.jl")
export MOGA_Real, MOGA_Binary, MOGA_Integer

include("core.jl")
export Problem, Algorithm, Solution, Archive

include("operators.jl")
export random_candidate, PM, SBX

include("algorithms.jl")
export NSGAII

end # module
