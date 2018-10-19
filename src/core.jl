using LinearAlgebra
import Base.copy

abstract type Algorithm end
abstract type GASolution end

#= Problem struct =#
"""
    Problem

Expects `eval_function` to have the form `f(s.x)` where `s` is a `Solution` and
to return a `Vector{Float64}`. If `has_constraint` is true, expects the final
value of the returned objective vector to be zero if solution is feasible, or
a positive Float64 indicating size of the constraint violation if infeasible.
"""
struct Problem
    nobjs::Int64
    maximize_objectives::Vector{Bool} # false -> minimize, true -> maximize
    nvars::Int64
    var_types::Vector{MOGA_Type}
    eval_fn::Function
    has_constraint::Bool
end

"no constraints constructor"
function Problem(nobjs::Int64, maximize_objectives::Vector{Bool},
                 nvars::Int64, var_types::Vector{MOGA_Type},
                 eval_fn::Function)
    Problem(nobjs, maximize_objectives, nvars, var_types,
            eval_fn, false)
end

#= Solutions =#
"""
    Solution

Type for candidate solutions, currently designed for NSGAII and NSGAIII.
"""
mutable struct Solution <: GASolution
    problem::Problem
    x::Vector{Float64}
    evaluated::Bool
    objectives::Vector{Float64}
    crowding_distance::Float64
    rank::Int64
    feasible::Bool
    constraint_violation::Float64
end

"no constraints constructor"
function Solution(problem::Problem,
                  x::Vector{Float64},
                  evaluated::Bool,
                  objectives::Vector{Float64},
                  crowding_distance::Float64,
                  rank::Int64)
    Solution(problem, x, evaluated, objectives, crowding_distance, rank,
             true, 0.0)
end

"allocates a new solution that copies argument"
function copy(s::Solution)
    return Solution(
        s.problem,
        copy(s.x),
        s.evaluated,
        copy(s.objectives),
        s.crowding_distance,
        s.rank,
        s.feasible,
        s.constraint_violation
    )
end

"Copies second argument solution into first argument"
function copyinto!(s_to::Solution, s_from::Solution)
    s_to.problem = s_from.problem
    s_to.x .= s_from.x
    s_to.evaluated = s_from.evaluated
    s_to.objectives .= s_from.objectives
    s_to.crowding_distance = s_from.crowding_distance
    s_to.rank = s_from.rank
    s_to.feasible = s_from.feasible
    s_to.constraint_violation = s_from.constraint_violation
end

"""
    evaluate!(s::Solution)

Calls `s.problem.eval_fn`, saves into `s.objectives`, sets `s.evaluated` to
true. If `s.problem.constraint_function` is defined, will also call this.
"""
function evaluate!(s::Solution)
    if !s.evaluated
        result = s.problem.eval_fn(s.x)
        if s.problem.has_constraint
            println("have constraint")
            s.feasible = result[end] == 0.0
            s.constraint_violation = result[end]
            s.objectives = result[1:end-1]
            @show(s.objectives)
        else
            println("no constraint")
            s.objectives = result
            @show(s.objectives)
        end
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

        if prob.maximize_objectives[i] # == MAXIMIZE
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


"""
    nondominated(solutions::Vector{Solution})

Returns all nondominated elements of `solutions`. This requires constructing
an `Archive` and inserting all solutions, so I'm not sure it's particularly fast.
"""
function nondominated(solutions::Vector{Solution})
    archive = Archive()
    insert_solutions!(archive, solutions)
    return archive.solutions
end


"""
    nondominated_sort(solutions::Vector{Solution})

Implements the nondominated sort for NSGAII. This function assigns front rank
and crowding distance scores to each solution in `solutions`. Note: better
function name might be `assign_rank_and_crowding!`, since it does not return
a sorted list.
"""
function nondominated_sort(solutions::Vector{Solution})
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


"""
    nondominated_cmp(x::Solution, y::Solution)

Returns -1 if x dominates y, 1 if y dominates x, and 0 otherwise. In the case
that both x and y are feasible, dominance is determined first by front rank,
and then by crowding distance. If only one solution is feasible, it automatically
dominates. If both are infeasbile, the solution with the smaller constraint
violation dominates.
"""
function nondominated_cmp(x::Solution, y::Solution)
    if x.feasible && y.feasible
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
    elseif x.feasible && !y.feasible
        return -1
    elseif !x.feasible && y.feasible
        return 1
    else
        if x.constraint_violation < y.constraint_violation
            return -1
        elseif y.constraint_violation < x.constraint_violation
            return 1
        else
            return 0
        end
    end
end


"""
    crowding_distance(solutions::Vector{Solution})

Assigns crowding distance, as defined in Deb 2002, to each element of
`solutions`.
"""
function crowding_distance(solutions::Vector{Solution})
    for s in solutions
        s.crowding_distance = 0.0
    end
    
    println("crowding_distance")
    @show(solutions[1].objectives)
    nobjs = length(solutions[1].objectives)
    @show(nobjs)

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


"""
    nondominated_truncate(pop::Vector{Solution}, num::Int64,
                          dominance::Function=nondominated_cmp)

Implements nondominated truncate used in NSGAII. `num` gives the desired number
of individuals in the truncated population. The `dominance` function is used to
determine sort order, so dominance(x, y) should return -1 if x dominates, 1 if
y dominates, and 0 if neither dominates.
"""
function nondominated_truncate(pop::Vector{Solution}, num::Int64,
                               dominance::Function=nondominated_cmp)
    # maybe more efficient to not sort the whole thing but pull out by
    # rank until reach the marginal rank class. Not sure. This is easy.
    k(c1, c2) = dominance(c1, c2) < 0
    sorted = sort(pop, lt=k)    # TODO: sort in place (sort!)?
    return sorted[1:num]
end


"""
    reference_point_truncate(
        pop::Vector{Solution}, N::Int64, refpts::Array{Float64,2})

Implements reference point truncation algorithm used in NSGAIII. `N` is the
desired population size, and `refpts` is an array containing the reference
points.
"""
function reference_point_truncate(pop::Vector{Solution}, N::Int64,
              refpts::Array{Float64, 2})

    M = length(pop[1].objectives)

    # find cutoff dominance rank
    ranks = [p.rank for p in pop]
    rankcounts = zeros(Int64, maximum(ranks))
    for r in ranks
        rankcounts[r] += 1
    end
    marginalrank = findfirst(x -> x >= N, cumsum(rankcounts))
    @show marginalrank

    St = [i for i=1:length(pop) if pop[i].rank <= marginalrank]

    if length(St) == N
        return [copy(pop[i]) for i in St]
    end
    # else:
    Pnext = [i for i in St if pop[i].rank < marginalrank]
    marginalfront = [i for i in St if pop[i].rank == marginalrank]

    ## Find ideal point: minimal value for each objective across all solutions
    zmin = zeros(M)
    for j = 1:M
        zmin[j] = minimum([pop[i].objectives[j] for i in St])
    end
    @show zmin

    ## FIND EXTREME VALUES
    ## the goal of this operation is to find, for each axis j, the solution that
    ## minimizes the value of non-j values. This isn't exactly the same as
    ## maximizing j if the solution space has a weird shape.
    ## TODO: I think this is right, but I'm not clear about the max vs mins here...
    expts = [-1 for j = 1:M]
    for j = 1:M  # looping through objectives/axes
        w = 0.00001 * ones(M)
        w[j] = 1.0
        maxval = -Inf
        for i in St
            val = maximum((pop[i].objectives) ./ w)
            if val > maxval
                expts[j] = i
                maxval = val
            end
        end
    end

    ## FIND INTERCEPT VALUES
    ## ENTER: Have found one extreme point per objective. These define a hyper-plane.
    ## EXIT: Find the intersections of the hyper-plane with each objective axis.
    ##    These will be used to rescale the objective scores for these axes.
    ##
    b = ones(M)
    Z = transpose(reduce(hcat, [pop[i].objectives for i in expts]))
    # try Z = factorize(Z)
    if rank(Z) == M
        a = 1.0 ./ (Z\b)
    else
        maxvals = zeros(M)
        for j=1:M
            maxvals[j] = maximum(pop[i].objectives[j] for i in St)
        end
        a = 1.0 ./ maxvals
    end

    # fp = [p.objectives - zmin for p in St] # ::Array{Array{Float64, 1}, 1}
    # fn = [f ./ a for f in fp] # ::Array{Array{Float64, 1}, 1}

    @show length(Pnext)
    @show length(marginalfront)

    fp = Array{Float64,1}[(pop[i].objectives - zmin) ./ a for i in Pnext]
    fm = Array{Float64,1}[(pop[i].objectives - zmin) ./ a for i in marginalfront]

    @show typeof(fp)
    @show typeof(fm)

    ## ASSOCIATE SOLUTIONS WITH REFERENCE POINTS
    ## Do some linear algebra to find distances between the objective values
    ## and lines passing through each reference point.
    fp_ref_index, fp_ref_distance = associate_points(fp, refpts)
    fm_ref_index, fm_ref_distance = associate_points(fm, refpts)
    fm_index_dict = Dict(i => d for (i, d) in zip(marginalfront, fm_ref_index))
    fm_dist_dict = Dict(i => d for (i, d) in zip(marginalfront, fm_ref_distance))


    ## NICHING PROCEDURE
    ## select points from the marginal front to include in P_{t+1}
    ## Prefer points associated with reference points with fewest associated
    ## solutions.

    K = N - length(Pnext)
    k = 1
    refpt_deck = Set(1:length(refpts))
    candidate_deck = Set(marginalfront) # set of indices into pop

    ρ = Dict(rp => 0 for rp in refpt_deck)
    for i in fp_ref_index
        ρ[i] += 1
    end

    while k <= K
        # find all reference points with minimal # of associated solutions from Pnext

        ρmin = minimum(values(ρ))
        # @show(ρmin)
        Jmin = [r for r in refpt_deck if ρ[r] == ρmin]
        # @show(Jmin)  # set of ref pts with minimal representation
        jbar = rand(Jmin)
        # find the set of all solutions in the marginal front associated with this ref point
        Ij = [i for i in candidate_deck if fm_index_dict[i] == jbar]
        if length(Ij) > 0
            # pick a point
            if ρ[jbar] == 0
                # pick the closest point
                distances = [fm_dist_dict[i] for i in Ij]
                chosen = Ij[argmin(distances)]
            else
                # choose any of the points
                chosen = rand(Ij)
            end
            push!(Pnext, chosen)
            ρ[jbar] += 1
            pop!(candidate_deck, chosen)
            k += 1
        else
            # case with no solutions in marginal front associated with this
            # reference point. Remove it from the pool for this iteration.
            pop!(refpt_deck, jbar)
            pop!(ρ, jbar)
        end
    end

    newpop = [copy(pop[i]) for i in Pnext] ## TODO: copies necessary?
    return newpop
end
