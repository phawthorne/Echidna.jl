using Echidna
using Plots
using Test

gr()

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
ndivs = 20
# refpts = generate_regular_reference_points(2, ndivs)
algo = NSGAIII(zdt1_problem, ZDT1, ndivs, niters)

seedgenomes = [(n/10.0)*ones(nvars) for n in 1:9]
seedpop = [candidate_from_genome(zdt1_problem, g) for g in seedgenomes]

pop = init_pop(algo, seedpop=seedpop)

## test one generation

# create new indivs
N = algo.population_size
for i = 1:(N/2)
    p1 = tournament_selector(pop, 2, dominance=nondominated_cmp)
    p2 = tournament_selector(pop, 2, dominance=nondominated_cmp)
    c1, c2 = [PM(c) for c in SBX(p1, p2)]
    for c in [c1, c2]
        c.objectives = algo.eval_fn(c.x)
        c.evaluated = true
    end
    push!(pop, c1, c2)
end

refpts = algo.reference_points

M = length(pop[1].objectives)

# find cutoff dominance rank
ranks = [p.rank for p in pop]
rankcounts = zeros(Int64, maximum(ranks))
for r in ranks
    rankcounts[r] += 1
end
marginalrank = findfirst(x -> x >= N, cumsum(rankcounts))

St = [i for i=1:length(pop) if pop[i].rank <= marginalrank]

if length(St) == N
    return St
end
# else:
Pnext = [i for i in St if pop[i].rank < marginalrank]
marginalfront = [i for i in St if pop[i].rank == marginalrank]

## Find ideal point: minimal value for each objective across all solutions
zmin = zeros(M)
for j = 1:M
    zmin[j] = minimum([pop[i].objectives[j] for i in St])
end

expts = [-1 for j = 1:M]
for j = 1:M  # looping through objectives/axes
    w = 0.00001 * ones(M)
    w[j] = 1.0
    maxval = -Inf
    for i in St
        # val = maximum((pop[i].objectives) ./ w)
        val = sum(pop[i].objectives ./ w)
        if val > maxval
            expts[j] = i
            maxval = val
        end
    end
end

extreme_points = [pop[i].objectives for i in expts]
x = [p.objectives[1] for p in pop]
y = [p.objectives[2] for p in pop]
scatter(x, y)
exx = [p[1] for p in extreme_points]
exy = [p[2] for p in extreme_points]
scatter!(exx, exy)

## FIND INTERCEPT VALUES
b = ones(nobjs)
Z = transpose(reduce(hcat, [pop[i].objectives for i in expts]))
a = 1.0 ./ (Z\b)

ax = [a[1], 0.0]
ay = [0.0, a[2]]
plot!(ax, ay)


fp = [(pop[i].objectives - zmin) ./ a for i in Pnext]
fm = [(pop[i].objectives - zmin) ./ a for i in marginalfront]

fpx = [p[1] for p in fp]
fpy = [p[2] for p in fp]
fmx = [p[1] for p in fm]
fmy = [p[2] for p in fm]
scatter(fpx, fpy)
scatter!(fmx, fmy)



## ASSOCIATE SOLUTIONS WITH REFERENCE POINTS
## Do some linear algebra to find distances between the objective values
## and lines passing through each reference point.
fp_ref_index, fp_ref_distance = associate_points(fp, refpts)
fm_ref_index, fm_ref_distance = associate_points(fm, refpts)
fm_index_dict = Dict(i => d for (i, d) in zip(marginalfront, fm_ref_index))
fm_dist_dict = Dict(i => d for (i, d) in zip(marginalfront, fm_ref_distance))


ρ = zeros(Int64, length(refpts))
for i in fp_ref_index
    ρ[i] += 1
end

K = N - length(Pnext)
k = 1
refpt_deck = Set(1:length(refpts))
candidate_deck = Set(marginalfront) # set of indices into pop


while k <= K
    # find all reference points with minimal # of associated solutions from Pnext
    ρmin = minimum(ρ)
    Jmin = [r for r in refpt_deck if ρ[r] == ρmin]  # set of ref pts with minimal representation
    jbar = rand(Jmin)
    @show jbar
    # find the set of all solutions in the marginal front associated with this ref point
    Ij = [i for i in candidate_deck if fm_index_dict[i] == jbar]
    @show Ij
    if length(Ij) > 0
        # pick a point
        if ρ[jbar] == 0
            # pick the closest point
            distances = [fm_dist_dict[i] for i in Ij]
            chosen = Ij[argmin(distances)]
        else
            # choose any of the points
            chosen = rand(Ij)
        end
        @show chosen
        push!(Pnext, chosen)
        ρ[jbar] += 1
        pop!(candidate_deck, chosen)
        k += 1
    else
        # case with no solutions in marginal front associated with this
        # reference point. Remove it from the pool for this iteration.
        pop!(refpt_deck, jbar)
    end
end


# result = garun(algo, seedpop=seed_indivs)
# @test length(result) == pop_size
