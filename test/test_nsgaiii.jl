using Echidna
using Test

print("Testing NSGAIII on ZDT1", "\n")
nobjs = 3
nvars = 30
ndivs = 5
pop_size = 100
n_iters = 10
zdt1_problem = Problem(
    2, [false for i in 1:nobjs], nvars,
    [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)


config = NSGAIII(zdt1_problem, ZDT1, nobjs, ndivs)
result = run(config)
@test length(result) == pop_size
