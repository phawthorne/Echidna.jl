include("core.jl")
include("algorithms.jl")

prob = Problem(2, [false, false])
sols = Vector{Solution}()
N = 10
for i=0:N
    r = i*(pi/2)/N
    objs = [cos(r), sin(r)]
    push!(sols, Solution(prob, [0.0], objs, 0.0, 0))
end

crowding_distance(sols)

archive = Archive(compare_pareto_dominance, Vector{Solution}())
insert_solutions!(archive, sols)
