using Parameters
include("core.jl")


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
