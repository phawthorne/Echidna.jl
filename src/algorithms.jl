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
        population = nsgaii_iteration(algo, population)
    end
    return population
end

function nsgaii_initialization(algo::NSGAII)
    N = algo.population_size
    population = [random_candidate(algo.problem) for i in 1:N]
    for p in population
        p.objectives = algo.eval_fn(p)
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
            c.objectives = algo.eval_fn(c)
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
