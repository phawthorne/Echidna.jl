using Echidna
using Test


# TEST ZDT1
x = [.1 for i in 1:30]
@show obj = ZDT1(x)
@test obj[1] == 0.1
@test obj[2] == 1.464110105645933

nobjs = 2
nvars = 30
zdt1_problem = Problem(
    2,
    [false for i in 1:nobjs],
    nvars,
    [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)
sol = Solution(
    zdt1_problem,
    [.1 for i in 1:nvars],
    false,
    [0.0, 0.0],
    0.0,
    0
)
# @show sol.objectives = ZDT1(sol)
@show evaluate!(sol)
@test sol.objectives[1] == 0.1
@test sol.objectives[2] == 1.464110105645933
@test sol.evaluated


#### TEST WITH TOURNAMENT
pop_size = 100
pop = [random_candidate(zdt1_problem) for i in 1:pop_size]
for i in 1:pop_size
    evaluate!(pop[i])
end
