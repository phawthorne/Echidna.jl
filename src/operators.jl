using Random

"""
    random_candidate(problem::Problem)

Returns a randomized `Solution` based on specification in `problem.var_types`.
"""
function random_candidate(problem::Problem)
    return Solution(
        problem,
        [rand(t) for t in problem.var_types],
        false,
        zeros(problem.nobjs),
        0.0,
        0
    )
end

"""
    candidate_from_genome(problem::Problem, genome)

Returns a Solution object that packs genome as the solution's decision variable.
"""
function candidate_from_genome(problem::Problem, genome)
    return Solution(
        problem,
        genome,
        false,
        zeros(problem.nobjs),
        0.0,
        0
    )
end

function tournament_selector(population::Vector{Solution}, tournament_size::Int;
                             dominance::Function=compare_pareto_dominance)
    popsize = length(population)

    winner = rand(population)

    for i in 1:tournament_size-1
        challenger = rand(population)
        if dominance(winner, challenger) > 0
            winner = challenger
        end
    end

    return winner
end


function PM(parent::Solution, probability=1, distribution_index=20.0)
    child = copy(parent)
    problem = child.problem

    if typeof(probability) <: Int
        probability /= sum([typeof(t) == MOGA_Real for t in problem.var_types])
    end

    for i in 1:problem.nvars
        if typeof(problem.var_types[i]) == MOGA_Real
            if rand() < probability
                child.x[i] = pm_mutation(child.x[i],
                                         problem.var_types[1].min_value,
                                         problem.var_types[1].max_value,
                                         distribution_index)
                child.evaluated = false
            end
        end
    end
    return child
end


function PM!(s::Solution, probability=1, distribution_index=20.0)
    problem = s.problem
    if typeof(probability) <: Int
        probability /= sum([typeof(t) == MOGA_Real for t in problem.var_types])
    end

    for i in 1:problem.nvars
        if typeof(problem.var_types[i]) == MOGA_Real
            if rand() < probability
                s.x[i] = pm_mutation(s.x[i],
                                     problem.var_types[1].min_value,
                                     problem.var_types[1].max_value,
                                     distribution_index)
                s.evaluated = false
            end
        end
    end
end


function pm_mutation(x::Float64, lb::Float64, ub::Float64, di::Float64)
    u = rand()
    dx = ub - lb

    if u < 0.5
        bl = (x - lb) / dx
        b = 2*u + (1-2*u)*(1-bl)^(di+1)
        delta = b^(1/(di + 1)) - 1
    else
        bu = (ub-x)/dx
        b = 2*(1-u) + 2(u-0.5)*(1-bu)^(di+1)
        delta = 1 - b^(1/(di+1))
    end
    x = x + delta*dx
    x = clip(x, lb, ub)

    return x
end


function SBX(p1::Solution, p2::Solution, probability=1.0, distribution_index=15.0)
    c1 = copy(p1)
    c2 = copy(p2)

    if rand() <= probability
        problem = c1.problem
        nvars = problem.nvars
        for i = 1:nvars
            if typeof(problem.var_types[i]) == MOGA_Real
                if rand() <= 0.5
                    x1 = c1.x[i]
                    x2 = c2.x[i]
                    lb = problem.var_types[i].min_value
                    ub = problem.var_types[i].max_value

                    x1, x2 = sbx_crossover(x1, x2, lb, ub, distribution_index)

                    c1.x[i] = x1
                    c2.x[i] = x2
                    c1.evaluated = false
                    c2.evaluated = false
                end
            end
        end
    end
    return c1, c2
end

"Note that this function modifies s1 and s2"
function SBX!(s1::Solution, s2::Solution, probability=1.0, distribution_index=15.0)
    if rand() <= probability
        problem = s1.problem
        nvars = problem.nvars
        for i = 1:nvars
            if typeof(problem.var_types[i]) == MOGA_Real
                if rand() <= 0.5
                    x1 = s1.x[i]
                    x2 = s2.x[i]
                    lb = problem.var_types[i].min_value
                    ub = problem.var_types[i].max_value

                    x1, x2 = sbx_crossover(x1, x2, lb, ub, distribution_index)

                    s1.x[i] = x1
                    s2.x[i] = x2
                    s1.evaluated = false
                    s2.evaluated = false
                end
            end
        end
    end
end

function sbx_crossover(x1::Float64, x2::Float64, lb::Float64, ub::Float64, di::Float64)
    dx = abs(x2 - x1)

    if dx > eps()
        if x2 > x1
            y1 = x1
            y2 = x2
        else
            y1 = x2
            y2 = x1
        end

        rx = rand()

        # calc x1
        b = 1.0 / (1.0 + (2.0 * (y1 - lb) / (y2 - y1)))
        alpha = 2.0 - b ^ (di + 1.0)
        if rx <= 1.0 / alpha
            alpha *= rx
            betaq = alpha ^ (1.0/(di+1.0))
        else
            alpha = 1.0 / (2.0 - alpha*rx)
            betaq = alpha ^ (1.0/(di+1.0))
        end
        x1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))

        # calc x2
        b = 1.0 / (1.0 + (2.0 * (ub - y2) / (y2 - y1)))
        alpha = 2.0 - b^(di + 1.0)
        if rx <= 1.0 / alpha
            alpha = alpha * rx
            betaq = alpha ^ (1.0 / (di + 1.0))
        else
            alpha = alpha * rx
            alpha = 1.0 / (2.0 - alpha)
            betaq = alpha ^ (1.0 / (di + 1.0))
        end
        x2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))

        # randomly swap the values
        if Random.bitrand()[1]
            x1, x2 = x2, x1
        end

        x1 = clip(x1, lb, ub)
        x2 = clip(x2, lb, ub)

    end
    return x1, x2
end

function clip(x::Real, lb::Real, ub::Real)
    return max(lb, min(x, ub))
end
