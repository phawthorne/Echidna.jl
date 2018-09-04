using LinearAlgebra

function find_extreme_points(points::Array{Vector{T}, 1},
                             zmin::Union{Missing, Vector{T}}=missing) where T <: Real
    M = length(points[1])

    @show(zmin)
    if ismissing(zmin)
        zmin = zeros(M)
    end
    @show zmin

    extremepts = [-1 for j = 1:M]
    for j = 1:M
        w = 0.00001 * ones(M)
        w[j] = 1.0
        minval = -Inf
        for (i, p) in enumerate(points)
            val = maximum((p-zmin) ./ w)
            if val > minval
                extremepts[j] = i
                minval = val
            end
        end
    end

    return extremepts
end


function find_hyperplane_params(points::Array{Vector{T}, 1}, b=1) where T <: Real
    M = length(points[1])

    interceptvals = zeros(M)

    Z = transpose(reduce(hcat, points))
    bv = b * ones(M)

    a = Z\bv
    return a
end

function find_axis_intercepts(points::Array{Vector{T}, 1}, b=1) where T <: Real
    return 1 ./ find_hyperplane_params(points, b)
end



function associate_points(
        canpoints::Array{Vector{T}, 1},
        refpoints::Array{Vector{T}, 1}) where T <: Real
    N = length(canpoints)
    R = length(refpoints)
    # for each candidate point, want to find the closest reference point,
    # measured as perpendicular distance to line through ref point and origin

    closest_ref_index = zeros(Int64, N)
    closest_ref_dist = zeros(Float64, N)
    dists = zeros(R)
    for n = 1:N
        s = canpoints[n]
        for r = 1:R
            w = refpoints[r]
            dists[r] = norm( s - (dot(w, s)/dot(w, w)) * w )
        end
        closest_ref_index[n] = argmin(dists)
        closest_ref_dist[n] = dists[closest_ref_index[n]]
    end
    return closest_ref_index, closest_ref_dist
end


function associate_points(
        canpoints::Array{Vector{T}, 1},
        refpoints::Array{T, 2}) where T <: Real
    N = length(canpoints)
    R = size(refpoints, 1)
    # for each candidate point, want to find the closest reference point,
    # measured as perpendicular distance to line through ref point and origin

    closest_ref_index = zeros(Int64, N)
    closest_ref_dist = zeros(Float64, N)
    dists = zeros(R)
    for n = 1:N
        s = canpoints[n]
        for r = 1:R
            w = refpoints[r, :]
            dists[r] = norm( s - (dot(w, s)/dot(w, w)) * w )
        end
        closest_ref_index[n] = argmin(dists)
        closest_ref_dist[n] = dists[closest_ref_index[n]]
    end
    return closest_ref_index, closest_ref_dist
end
