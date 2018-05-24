import Base.copy
include("var_types.jl")


struct Problem
    nobjs::Int64
    directions::Vector{Bool}
    nvars::Int64
    var_types::Vector{MOGA_Type}
end

abstract type Algorithm end

# Solutions
mutable struct Solution
    problem::Problem
    x::Vector{Float64}
    evaluated::Bool
    objectives::Vector{Float64}
    crowding_distance::Float64
    rank::Int64
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

# Archives
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


# Helper functions
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


function nondominated(solutions)
    archive = Archive()
    insert_solutions!(archive, solutions)
    return archive.solutions
end

function nondominated_sort(solutions)
    rank::Int64 = 0

    for s in solutions
        s.rank = -1
    end

    while length(solutions) > 0
        archive = Archive(compare_pareto_dominance, Vector{Solution}())
        insert_solutions!(archive, solutions)

        for s in archive.solutions
            s.rank = rank
        end

        crowding_distance(archive.solutions)
        solutions = [s for s in solutions if s.rank == -1]
        rank += 1
    end
end

function nondominated_cmp(x::Solution, y::Solution)
    if x.rank == y.rank
        if -x.crowding_distance < -y.crowding_distance
            return -1
        elseif -x.crowding_distance > -y.crowding_distance
            return 1
        else
            return 0
        end
    else
        if x.rank < y.rank
            return -1
        elseif x.rank > y.rank
            return 1
        else
            return 0
        end
    end
end

function crowding_distance(solutions)
    for s in solutions
        s.crowding_distance = 0.0
    end

    nobjs = length(solutions[1].objectives)

    if length(solutions) < 3
        for s in solutions
            s.crowding_distance[1] = Inf
        end
    else
        # main case
        for i in 1:nobjs
            sorted_solutions = sort!(solutions, by=s->s.objectives[i])
            min_value = sorted_solutions[1].objectives[i]
            max_value = sorted_solutions[end].objectives[i]

            sorted_solutions[1].crowding_distance += Inf
            sorted_solutions[end].crowding_distance += Inf

            for j in 2:length(solutions)-1
                if max_value - min_value < eps(Float64)
                    sorted_solutions[j].crowding_distance = Inf
                else
                    diff = sorted_solutions[j+1].objectives[i] -
                           sorted_solutions[j-1].objectives[i]
                    sorted_solutions[j].crowding_distance += diff /
                                                (max_value - min_value)
                end
            end
        end
    end
end
