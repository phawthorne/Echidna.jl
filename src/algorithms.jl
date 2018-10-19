using Parameters

@with_kw mutable struct NSGAII <: Algorithm
    problem::Problem
    eval_fn::Function
    population_size::Int64
    n_iters::Int64
    archive::Archive
    archive_frequency::Int64
end

function NSGAII(problem::Problem, eval_fn::Function, pop_size::Int64, n_iters::Int64)
    archive = Archive(compare_pareto_dominance, Vector{Solution}())
    archive_frequency = 100000000
    return NSGAII(problem, eval_fn, pop_size, n_iters, archive, archive_frequency)
end

function garun(algo::NSGAII; seedpop::Vector{Solution}=Vector{Solution}())
    population = init_pop(algo; seedpop=seedpop)

    for i in 1:algo.n_iters
        if i % algo.archive_frequency == 0
            insert_solutions!(algo.archive, population)
        end
        # println("$i")
        population = iter_generation(algo, population)
    end
    return population
end

function init_pop(algo::NSGAII; seedpop::Vector{Solution}=Vector{Solution}())
    N = algo.population_size
    population = vcat(seedpop,
                      [random_candidate(algo.problem) for i in 1:(N-length(seedpop))])

    for p in population
        evaluate!(p)
    end

    nondominated_sort(population)
    return population
end

function iter_generation(algo::NSGAII, population::Vector{Solution})
    # create N new indivs: binary tournament -> SBX/PM
    N = algo.population_size
    for i = 1:(N/2)
        p1 = tournament_selector(population, 2, dominance=nondominated_cmp)
        p2 = tournament_selector(population, 2, dominance=nondominated_cmp)
        c1, c2 = [PM(c) for c in SBX(p1, p2)]
        for c in [c1, c2]
            evaluate!(c)
        end
        push!(population, c1, c2)
    end

    # assign ranks and crowding_distance
    nondominated_sort(population)

    # select best N as next generation
    population = nondominated_truncate(population, N)

    return population
end

# TODO: add Archive to NSGAIII
@with_kw mutable struct NSGAIII <: Algorithm
    problem::Problem
    eval_fn::Function
    nobjs::Int64
    ndivs::Int64
    population_size::Int64
    n_iters::Int64
    reference_points::Array{Float64, 2}
end

function NSGAIII(problem::Problem, evalfn::Function, ndivs::Int, niters::Int)
    refpts = generate_regular_reference_points(problem.nobjs, ndivs)
    return NSGAIII(
        problem,
        evalfn,
        problem.nobjs,
        ndivs,
        Int(ceil(size(refpts, 1)/4)*4),
        niters,
        refpts
    )
end

function garun(algo::NSGAIII; seedpop::Vector{Solution}=Vector{Solution}())
    println("Echidna.garun")

    population = init_pop(algo; seedpop=seedpop)

    for g in 1:algo.n_iters
        println("")
        println("generation: ", g)
        population = iter_generation(algo, population)
    end

    return population
end

function init_pop(algo::NSGAIII; seedpop::Vector{Solution}=Vector{Solution}())
    N = algo.population_size
    population = vcat(seedpop,
                      [random_candidate(algo.problem) for i in 1:(N-length(seedpop))])

    for p in population
        evaluate!(p)
    end

    nondominated_sort(population)

    return population
end

function iter_generation(algo::NSGAIII, population::Vector{Solution})
    # TODO: this is copied from NSGAII - make sure there shouldn't be other changes
    # create N new indivs: binary tournament -> SBX/PM
    N = algo.population_size
    for i = 1:(N/2)
        p1 = tournament_selector(population, 2, dominance=nondominated_cmp)
        p2 = tournament_selector(population, 2, dominance=nondominated_cmp)
        c1, c2 = [PM(c) for c in SBX(p1, p2)]
        for c in [c1, c2]
            evaluate!(c)
        end
        push!(population, c1, c2)
    end

    # assign ranks and crowding_distance
    nondominated_sort(population)

    # select best N as next generation
    population = reference_point_truncate(population, N, algo.reference_points)

    return population
end
