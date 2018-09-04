using Echidna
using Plots

gr()

f1(θ) = [cos(θ), sin(θ)]

pts = [f1(π/2 * x) for x in rand(10)]
x = [p[1] for p in pts]
y = [p[2] for p in pts]
scatter(x, y)

zmin = [1.0, 1.0]
extremepts = find_extreme_points(pts, zmin)

xx = [x[i] for i in extremepts]
yy = [y[i] for i in extremepts]
scatter!(xx, yy)


f2(x) = [x, exp(-x)]
pts = [f2(5*x) for x in rand(10)]
x = [p[1] for p in pts]
y = [p[2] for p in pts]
scatter(x, y)

extremepts = find_extreme_points(pts)

xx = [x[i] for i in extremepts]
yy = [y[i] for i in extremepts]
scatter!(xx, yy)

zpts = [pts[i] for i in extremepts]
find_axis_intercepts(zpts)



#### TESTING POINT ASSOCIATION
canpoints = [[0.5, 1.0], [2.0, 0.0]]
refpoints = [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]

index, dists = associate_points(canpoints, refpoints)
