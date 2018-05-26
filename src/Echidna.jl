module Echidna

include("var_types.jl")
export MOGA_Real, MOGA_Binary, MOGA_Integer

include("core.jl")
export Problem, Algorithm, Solution, Archive, evaluate!,
       insert_solutions!, insert_solution!, compare_pareto_dominance

include("operators.jl")
export random_candidate, PM, SBX

include("algorithms.jl")
export NSGAII, run_nsgaii

include("problems.jl")
export ZDT1

end # module
