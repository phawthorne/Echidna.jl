

# TODO: this interface is a bit annoying - I wonder if there's a way to have
# a single function call that will work here.
function ZDT1(sol::Solution)
    f1 = sol.x[1]
    g = 1 + (9/(length(sol.x)-1)) * sum(sol.x[2:end])
    h = 1-sqrt(f1/g)
    return [f1, g*h]
end

function ZDT1(x::Vector{Float64})
    f1 = x[1]
    g = 1 + (9/(length(x)-1)) * sum(x[2:end])
    h = 1-sqrt(f1/g)
    return [f1, g*h]
end
