include("core.jl")
include("var_types.jl")

function random_candidate(problem::Problem)
    return Solution(
        problem,
        [rand(t) for t in problem.var_types],
        false
        zeros(problem.nobjs),
        0.0,
        0
    )
end

function PM(parent, probability=1, distribution_index=20.0)
    child = copy(parent)
    problem = child.problem

    if typeof(probability) <: Int
        probability /= sum([typeof(t) == MOGA_Real for t in problem.var_types])
    end

    for i in range(problem.nvars)
        if typeof(problem.var_types[i]) == MOGA_Real
            if rand() < probability
                child.x[i]
            end
        end
    end

end

function pm_mutation(x, lb, ub)

end
