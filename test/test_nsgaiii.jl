using Echidna
# using Plots
using Test

# gr()

println("Setting up ZDT1 test problem")
nobjs = 2
nvars = 30
zdt1_problem = Problem(
    nobjs, [false for i in 1:nobjs], nvars,
    [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)

println("Setting up NSGAIII config")
niters = 200
ndivs = 25
# refpts = generate_regular_reference_points(2, ndivs)
algo = NSGAIII(zdt1_problem, ZDT1, ndivs, niters)
result = garun(algo)

@show result


# rx = [s.objectives[1] for s in result]
# ry = [s.objectives[2] for s in result]
# scatter(rx, ry)

# seedgenomes = [(n/10.0)*ones(nvars) for n in 1:9]
# seedpop = [candidate_from_genome(zdt1_problem, g) for g in seedgenomes]
