include("core.jl")
include("var_types.jl")

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

function PM(parent, probability=1, distribution_index=20.0)
    child = copy(parent)
    problem = child.problem

    if typeof(probability) <: Int
        probability /= sum([typeof(t) == MOGA_Real for t in problem.var_types])
    end

    for i in range(problem.nvars)
        if typeof(problem.var_types[i]) == MOGA_Real
            if rand() < probability
                child.x[i] = pm_mutation(child.x[i],
                                         problem.var_types[].min_value,
                                         problem.var_types[].max_value,
                                         distribution_index)
                child.evaluated = false
            end
        end
    end
    return child
end

function pm_mutation(x, lb, ub, di)
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

function SBX(p1, p2, probability=1.0, distribution_index=15.0)
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

                    x1, x2 = sbx_crossover(x1, x2, lu, ub, distribution_index)

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

function sbx_crossover(x1, x2, lb, ub, di)
    dx = x2 - x1

    if dx > eps()
        if x2 > x1
            y1 = x1
            y2 = x2
        else
            y1 = x2
            y2 = x1
        end

        beta = 1.0 / (1.0 + (2.0 * (y1 - lb) / (y2 - y1)))
        alpha = 2.0 - beta ^ (di + 1.0)

        rx = rand()

        if rx <= 1.0 / alpha
            alpha *= rx
            betaq = alpha ^ (1.0/(di+1.0))
        else
            alpha = 1.0 / (2.0 - alpha*rx)
            betaq = alpha ^ (1.0/(di+1.0))
        end

        x1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
        beta = 1.0 / (1.0 + (2.0 * (ub - y2) / (y2 - y1)))
        alpha = 2.0 - beta^(di + 1.0)

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
        if bitrand()[1]
            x1, x2 = x2, x1
        end

        x1 = clip(x1, lb, ub)
        x2 = clip(x2, lb, ub)

        return x1, x2
    end
end

function clip(x, lb, ub)
    return max(lb, min(x, ub))
end
