using Parameters

@with_kw mutable struct NSGAII <: Algorithm
    problem::Problem
    eval_fn::Function
    population_size::Int64
    n_iters::Int64
end

function run_nsgaii(algo::NSGAII)
    population = nsgaii_initialization(algo)

    for i in 1:algo.n_iters
        print("$i\n")
        population = nsgaii_iteration(algo, population)
    end
    return population
end

function nsgaii_initialization(algo::NSGAII)
    N = algo.population_size
    population = [random_candidate(algo.problem) for i in 1:N]
    population[1].x .= 0.0001   # This is a hack for the NNM GA
                                # TODO: allow some seeding

    for p in population
        p.objectives = algo.eval_fn(p.x)
        p.evaluated = true
    end
    # assign ranks and crowding_distance
    nondominated_sort(population)
    return population
end

function nsgaii_iteration(algo::NSGAII, population::Vector{Solution})
    # create N new indivs: binary tournament -> SBX/PM
    N = algo.population_size
    for i = 1:(N/2)
        p1 = tournament_selector(population, 2, dominance=nondominated_cmp)
        p2 = tournament_selector(population, 2, dominance=nondominated_cmp)
        c1, c2 = [PM(c) for c in SBX(p1, p2)]
        for c in [c1, c2]
            c.objectives = algo.eval_fn(c.x)
            c.evaluated = true
        end
        push!(population, c1, c2)
    end

    # assign ranks and crowding_distance
    nondominated_sort(population)

    # select best N as next generation
    population = nondominated_truncate(population, N)

    return population
end

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

function run(algo::NSGAIII)
    nobjs = 3
    ndivs = 10
    algo.reference_points = generate_regular_reference_points(nobjs, ndivs)
    nrps = size(algo.reference_points)[1]
    algo.population_size = Int(ceil(nrps/4)*4)

    population = init_pop(algo)

    for g in 1:algo.n_iters
        println("generation: ", g)
        population = iter_generation(algo, population)
    end

    return population
end

function init_pop(algo::NSGAIII; seedpop::Vector{GASolution}=Vector{GASolution}())
    N = algo.population_size
    population = vcat(seedpop,
                      [random_candidate(algo.problem) for i in 1:(N-length(seedpop))])

    for p in population
        p.objectives = algo.eval_fn(p.x)
        p.evaluated = true
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
            c.objectives = algo.eval_fn(c.x)
            c.evaluated = true
        end
        push!(population, c1, c2)
    end

    # assign ranks and crowding_distance
    nondominated_sort(population)

    # select best N as next generation
    population = reference_point_truncate(population, N)

    return population
end
