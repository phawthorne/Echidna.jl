using CSV
using DataFrames


function log_population(population::Vector{Solution}, gen, logdest)
    popsize = length(population)
    nobjs = population[1].problem.nobjs

    
