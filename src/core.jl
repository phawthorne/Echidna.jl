import Base.copy

struct Problem
    nobjs::Int64
    directions::Vector{Bool}    # false -> minimize, true -> maximize
    nvars::Int64
    var_types::Vector{MOGA_Type}
    eval_fn::Function
end

abstract type Algorithm end
abstract type GASolution end

# Solutions
mutable struct Solution <: GASolution
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

function evaluate!(s::Solution)
    if !s.evaluated
        s.objectives = s.problem.eval_fn(s.x)
        s.evaluated = true
    end
    return s.objectives
end

# Archives
mutable struct Archive
    dominance::Function
    solutions::Vector{Solution}
end

function insert_solution!(archive::Archive, solution::Solution)
    flags = [archive.dominance(solution, arch_sol) for
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
    rank::Int64 = 1  # TODO 08/30/18 - I just changed this to 1 to be Julian - hope nothing breaks

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
            s.crowding_distance = Inf
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

function nondominated_truncate(pop::Vector{Solution}, num::Int64,
                               dominance::Function=nondominated_cmp)
    # maybe more efficient to not sort the whole thing but pull out by
    # rank until reach the marginal rank class. Not sure. This is easy.
    k(c1, c2) = dominance(c1, c2) < 0
    sorted = sort(pop, lt=k)
    return sorted[1:num]
end

function reference_point_truncate(pop::Vector{Solution}, N::Int64,
                                  refpts::Array{Float64, 2})
    # find cutoff dominance rank
    ranks = [p.rank for p in pop]
    rankcounts = zeros(Int64, maximum(ranks))
    for r in ranks
        rankcounts[r] += 1
    end
    marginalrank = findfirst(x -> x >= N, cumsum(rankcounts))

    St = [p for p in pop if p.rank <= marginalrank]
    topick = N - length(St)
    if topick == 0
        return St
    end
    # else:
    marginalfront = [p for p in newpop if p.rank == marginalrank]

    ## Get normalized values
    M = length(St[1].objectives)
    zmin = zeros(M)
    for j = 1:M
        zmin = minimum([p.objectives[j] for p in St])
    end
    fp = [p.objectives - zmin for p in St]

    ## FIND EXTREME VALUES
    ## the goal of this operation is to find, for each axis j, the solution that
    ## minimizes the value of non-j values. This isn't exactly the same as
    ## maximizing j if the solution space has a weird shape.
    ##
    ## EXIT: extremepts contains the indices of solutions in St that are the
    ## extreme points for each objective axis.
    ##
    ## TODO: I'm not clear about the max vs mins here...
    extremepts = [-1 for j = 1:M]
    for j = 1:M
        w = 0.00001 * ones(M)
        w[j] = 1.0
        minval = -Inf
        for (i, p) in enumerate(St)
            val = maximum((p.objectives - zmin) ./ w)
            if val > minval
                extremepts[j] = i
                minval = val
            end
        end
    end

    ## FIND INTERCEPT VALUES
    ## ENTER: Have found one extreme point per objective. These define a hyper-plane.
    ## EXIT: Find the intersections of the hyper-plane with each objective axis.
    ##    These will be used to rescale the objective scores for these axes.
    ##
    b = ones(nobjs)
    Z = zeros(nobjs, nobjs)
    for r in 1:nobjs
        index = extremepts[r]
        for c in 1:nobjs
            Z[r,c] = fp[r,c]
        end
    end
    a = Z\b

    fn = [f ./ a for f in fp]

    ## ASSOCIATE SOLUTIONS WITH REFERENCE POINTS
    ## Do some linear algebra to find distances between the objective values
    ## and lines passing through each reference point.


    ## NICHING PROCEDURE
    ## select points from the marginal front to include in P_{t+1}
    ## Prefer points associated with reference points with fewest associated
    ## solutions.


    return pop
end
