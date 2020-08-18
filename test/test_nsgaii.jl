using Echidna
using Test

print("Testing NSGAII on ZDT1", "\n")
nobjs = 2
nvars = 30
pop_size = 100
n_iters = 10
zdt1_problem = Problem(
    2, [false for i in 1:nobjs], nvars,
    [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)

config = NSGAII(zdt1_problem, ZDT1, pop_size, n_iters)
result = garun(config)
@test length(result) == pop_size
