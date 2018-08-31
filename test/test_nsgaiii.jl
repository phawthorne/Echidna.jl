using Echidna
using Test

println("Setting up ZDT1 test problem")
nobjs = 2
nvars = 30
pop_size = 100
niters = 10
zdt1_problem = Problem(
    nobjs, [false for i in 1:nobjs], nvars,
    [MOGA_Real(0.0, 1.0) for i in 1:30],
    ZDT1
)

println("Setting up NSGAIII config")
ndivs = 50
config = NSGAIII(zdt1_problem, ZDT1, ndivs, niters)

seed_genomes = [(n/10.0)*ones(nvars) for n in 1:9]
seed_indivs = [candidate_from_genome(zdt1_problem, g) for g in seed_genomes]

result = garun(config, seedpop=seed_indivs)
@test length(result) == pop_size
