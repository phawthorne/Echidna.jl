using Echidna
using Base.Test

# test for compare_pareto_dominance
nobjs = 2
nvars = 30
zdt1_problem = Problem(
    2, [false for i in 1:nobjs], nvars, [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)

objective_vals = [0.,3.], [1., 1.], [2., 2.], [3., 0.]
sols = [ Solution(
    zdt1_problem, [.1 for i in 1:nvars], false, ov, 0.0, 0
) for ov in objective_vals]

expectation = [[0,0,0,0] [0,0,-1,0] [0,1,0,0] [0,0,0,0]]
for i=1:4
    for j=1:4
        @test compare_pareto_dominance(sols[i], sols[j]) == expectation[j, i]
    end
end

a = Archive(compare_pareto_dominance, Vector{Solution}())
insert_solutions!(a, sols)
@test length(a.solutions) == 3

# SWITCH OBJECTIVE GOALS
zdt1_problem = Problem(
    2, [true, false], nvars, [MOGA_Real(0.0, 1.0) for i in 1:30], ZDT1
)
sols = [ Solution(
    zdt1_problem, [.1 for i in 1:nvars], false, ov, 0.0, 0
) for ov in objective_vals]

expectation = [[0,1,1,1] [-1,0,0,1] [-1,0,0,1] [-1,-1,-1,0]]
for i=1:4
    for j=1:4
        @test compare_pareto_dominance(sols[i], sols[j]) == expectation[j, i]
    end
end

a = Archive(compare_pareto_dominance, Vector{Solution}())
insert_solutions!(a, sols)
@test length(a.solutions) == 1
