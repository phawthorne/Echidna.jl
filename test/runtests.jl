using Echidna
using Test

@testset "Echidna" begin
    @testset "Test Problems" begin include("test_problems.jl"); end
    @testset "Core" begin include("test_core.jl"); end
    @testset "NSGA-II" begin include("test_nsgaii.jl"); end
end

# include("test_nsgaiii.jl")
