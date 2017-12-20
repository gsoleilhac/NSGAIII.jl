using NSGAIII
using vOptGeneric, GLPKMathProgInterface

using PyPlot

function plot_pop(P)
    clf()
    pop = filter(x -> x.CV ≈ 0, P)
    non_dom = NSGAIII.fast_non_dominated_sort!(P)[1]
    pop = setdiff(pop, non_dom)
    p = plot3D(map(x -> x.y[1], pop), map(x -> x.y[2], pop), map(x -> x.y[3], pop), "bo", markersize=1)
    p = plot3D(map(x -> x.y[1], non_dom), map(x -> x.y[2], non_dom), map(x -> x.y[3], non_dom), "go", markersize=1)
    ax = gca()
    ax[:set_xlim]([-95, 110])
    ax[:set_ylim]([-60, 45])
    ax[:set_zlim]([-50, 110])
    !isinteractive() && show()
    sleep(0.1)
end


m = vModel(solver = GLPKSolverMIP())

@variable(m, 0 <=x[1:5] <= 10)
@variable(m, δ[1:3], Bin)

@addobjective(m, Max, dot([17,-12,-12,-19,-6], x) + dot([-73, -99, -81], δ))
@addobjective(m, Max, dot([2,-6,0,-12,13], x) + dot([-61,-79,-53], δ))
@addobjective(m, Max, dot([-20,7,-16,0,-1], x) + dot([-72,-54,-79], δ))

@constraint(m, sum(δ) <= 1)
@constraint(m, -x[2] + 6x[5] + 25δ[1] <= 52)
@constraint(m, -x[1] + 18x[4] + 18x[5] + 8δ[2] <= 77)
@constraint(m, 7x[4] + 9x[5] + 19δ[3] <= 66)
@constraint(m, 16x[1] + 20x[5] <= 86)
@constraint(m, 13x[2] + 7x[4] <= 86)


print(m)


res = nsga(0, 50, 30, m, fplot = plot_pop, pmut=0.3);
solve(m, method=:lex)

sort(res, by = x -> x[2][1])[1]
sort(res, by = x -> x[2][2])[1]
sort(res, by = x -> x[2][3])[1]
getY_N(m)
