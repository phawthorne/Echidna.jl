using Base.Filesystem
using CSV
using DataFrames

"Save objective values of a population to a file"
function log_population(population::Vector{Solution}, gen, log_dir)
    if ~isdir(log_dir)
        mkpath(log_dir)
    end
    
    popsize = length(population)
    nobjs = population[1].problem.nobjs

    df = DataFrame()
    df[!, "iteration"] = [gen for p = 1:popsize]
    for i=1:nobjs
        colname = "obj_$i"
        df[!, colname] = [population[p].objectives[i] for p = 1:popsize]
    end
    df[!, "born"] = [population[p].generation for p = 1:popsize]
    df[!, "rank"] = [population[p].rank for p = 1:popsize]

    output_file = Base.Filesystem.joinpath(log_dir, "log_$gen.csv")
    CSV.write(output_file, df)
    @show(gen)

end
    
