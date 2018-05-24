using Parameters

@with_kw mutable struct NSGAII <: Algorithm
    problem::Problem
    eval_fn::Function
    population_size::Int64
    population::Vector{Solution}
    iteration::Int64

end
