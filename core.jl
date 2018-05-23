import Base.copy
include("var_types.jl")

struct Problem
    nobjs::Int64
    directions::Vector{Bool}
    nvars::Int64
    var_types::Vector{MOGA_Type}
end

mutable struct Solution
    problem::Problem
    x::Vector{Float64}
    evaluated::Bool
    objectives::Vector{Float64}
    crowding_distance::Float64
    rank::Int64
end

mutable struct Archive
    dominance::Function
    solutions::Vector{Solution}
end


function insert_solution!(archive::Archive, solution::Solution)
    flags = [archive.dominance(arch_sol, solution) for
             arch_sol in archive.solutions]
    dominates = [x > 0 for x in flags]
    nondominated = [x == 0 for x in flags]

    if any(dominates)
        return false
    else
        archive.solutions = archive.solutions[nondominated]
        push!(archive.solutions, solution)
        return true
    end
end

function insert_solutions!(archive::Archive, solutions::Vector{Solution})
    for s in solutions
        insert_solution!(archive, s)
    end
end


"""
    compare_pareto_dominance(s1, s2)

returns -1 if `s1` dominates `s2`, 1 if `s2` dominates `s1`, and 0 otherwise
"""
function compare_pareto_dominance(s1, s2)
    prob = s1.problem
    # THIS IS WHERE TO INSERT CONSTRAINT CHECK

    dom1 = false
    dom2 = false

    for i in 1:prob.nobjs
        o1 = s1.objectives[i]
        o2 = s2.objectives[i]

        if prob.directions[i] # == MAXIMIZE
            o1 = -1 * o1
            o2 = -1 * o2
        end

        if o1 < o2
            dom1 = true
            if dom2
                return 0
            end
        elseif o1 > o2
            dom2 = true
            if dom1
                return 0
            end
        end
    end

    if dom1 == dom2
        return 0
    elseif dom1
        return -1
    else
        return 1
    end

end

function copy(s::Solution)
    return Solution(
        s.problem,
        copy(s.x),
        s.evaluated,
        copy(s.objectives),
        s.crowding_distance,
        s.rank
    )
end
