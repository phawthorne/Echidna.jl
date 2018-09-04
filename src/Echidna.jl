module Echidna

include("var_types.jl")
export MOGA_Real, MOGA_Binary, MOGA_Integer

include("core.jl")
export Problem, Algorithm, Solution, Archive, evaluate!,
       insert_solutions!, insert_solution!, compare_pareto_dominance,
       nondominated_cmp, nondominated_sort, crowding_distance,
       nondominated_truncate, reference_point_truncate

include("operators.jl")
export random_candidate, candidate_from_genome, PM, SBX, tournament_selector

include("algorithms.jl")
export NSGAII, run_nsgaii, NSGAIII, garun, init_pop, iter_generation

include("problems.jl")
export ZDT1

include("weights.jl")
export generate_regular_reference_points, generate_recursive

include("geom_ops.jl")
export find_extreme_points, find_hyperplane_params, find_axis_intercepts,
       associate_points

end # module
