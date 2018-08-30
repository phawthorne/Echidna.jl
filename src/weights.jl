

# def normal_boundary_weights(nobjs, divisions_outer, divisions_inner=0):
#     def generate_recursive(weights, weight, left, total, index):
#         if index == nobjs - 1:
#             weight[index] = float(left) / float(total)
#             weights.append(copy.copy(weight))
#         else:
#             for i in range(left+1):
#                 weight[index] = float(i) / float(total)
#                 generate_recursive(weights, weight, left-i, total, index+1)
#
#     def generate_weights(divisions):
#         weights = []
#         generate_recursive(weights, [0.0]*nobjs, divisions, divisions, 0)
#         return weights
#
#     weights = generate_weights(divisions_outer)
#
#     if divisions_inner > 0:
#         inner_weights = generate_weights(divisions_inner)
#
#         for i in range(len(inner_weights)):
#             weight = inner_weights[i]
#
#             for j in range(len(weight)):
#                 weight[j] = (1.0 / nobjs + weight[j]) / 2.0
#
#             weights.append(weight)
#
#     return weights


"""
    generate_regular_reference_points(M, p)

Generates reference points following the Das and Dennis systematic approach.
Notation follows Deb and Jain (NSGA-III paper)
M = number of objectives
p = number of divisions
"""
function generate_regular_reference_points(M, p)
    weights = zeros(0,M)
    weights = generate_recursive(weights, zeros(1, M), p, p, M, 1)
    return weights
end


function generate_recursive(weights, weight, left, total, M, index)
    if index == M
        # termination case
        weight[index] = left / total
        weights = [weights; weight]
    else
        for i in 0:left
            weight[index] = i / total
            weights = generate_recursive(weights, weight, left-i, total, M, index+1)
        end
    end
    return weights
end
