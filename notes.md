# Notes on reducing allocation of new pops

## Functions that allocate Solutions:
* algorithms.init_pop
    - calls operators.random_candidate
* algorithms.iter_generation
    - calls operators.PM and operators.SBX
* operators.random_candidate
* operators.candidate_from_genome

## Functions that copy:
* operators.PM
* operators.SBX
* core.reference_point_truncate

## Functions that discard:
* core.nondominated_truncate
* core.reference_point_truncate

## Implementation plan
* Need to add a Nursery [DataStructures.Circular.Deque] allocation to Algorithm structure and constructor.
* Right now random_candidate and candidate_from_genome only take a Problem, not an Algorithm, and they construct a new Solution. Should they require an Algorithm that contains a Nursery? I could also have a versions that take a Nursery. Need a way to copy an allocated Solution into the Nursery so that we don't accidentally explode the number of Solutions.
* The truncate functions need to be amended to stick non-selected Solutions back into the Nursery. Pretty straightforward for nondominated_truncate (i.e. sorted[num+1:end]). reference_point_truncate will require some thinking.
* Could rewrite the \_truncate functions to return list of indices for included pops. Then could just iterate through indices and either add to next pop or "decompose".

# Questions for self
* I think I'm checking for var_type here because PM only applies to real-valued genome elements. If we had a mixed type genome, would we sequentially apply PM and then another?


# Other notes:
* Probably don't want to deepcopy in operators.PM - this is creating copies of the problem, which includes a lot of things that don't need to be duplicated. Why didn't I use the custom copy(s::Solution)?
* No structure is actually holding the population right now. It gets created in garun() and acted on by iter_generation. Maybe I want to add it to the
* Merge cases where code is identical between algorithms? Use Union types.
* Add mutator and crossover functions to algorithm structure with PM and SBX as defaults.
* Can merge NSGAII and NSGAIII iter_generations by adding truncator function to the algorithm def (NSGAII: nondominated_truncate, NSGAIII: reference_point_truncate)
* Journal article Ahmadi 2015 uses hypervolume (Lebesgue measure) with a reference point as "goodness" measure of a frontier. Should implement this.


# Perf results
baseline for SeedMix:
`@btime results = garun(algo, seedpop=seedpop)
  1.668 s (10893234 allocations: 604.76 MiB)`

swap deepcopy for copy:
`@btime results = garun(algo, seedpop=seedpop);
  1.251 s (10729531 allocations: 559.29 MiB)`


```julia
    s1 = results[1]
    @btime s2 = copy(s1)
        328.455 ns (3 allocations: 3.94 KiB)
    @btime s2 = deepcopy(s1)
        32.595 μs (14 allocations: 8.34 KiB)
    @btime copyinto!(s2, s1)
        51.663 ns (0 allocations: 0 bytes)
```
