# NSGAII

[![Build Status](https://travis-ci.org/gsoleilhac/NSGAII.jl.svg?branch=master)](https://travis-ci.org/gsoleilhac/NSGAIII.jl)

[![codecov.io](http://codecov.io/github/gsoleilhac/NSGAII.jl/coverage.svg?branch=master)](http://codecov.io/github/gsoleilhac/NSGAIII.jl?branch=master)

## Installation
`Pkg.clone("https://github.com/gsoleilhac/NSGAIII.jl")`


## Usage : 

```
popsize = 200
nbGenerations = 100

#Define the number of division along each objective to generate the reference directions.
H = 5
#Alternatively, you can directly pass the reference directions as a Vector{Vector{Float64}} : 
#With two objectives, H = 2 is equivalent to 
H = [[1., 0.], [0.5, 0.5], [0., 1.]]

#define how to generate a random genotype : 
init_function = () -> randperm(N) #e.g : a permutation coding

#define how to evaluate a genotype : 
z(x) = z1(x), z2(x) ... # Must return a Tuple

nsga(popsize, nbGenerations, init_function, z, H)

#A constraint violation function can be passed with the keyword fCV
#It should return 0. if the solution is feasible and a value > 0 otherwise.

#Mutation probability can be changed with the keyword pMut (default is 5%)

nsga(popsize, nbGen, init_fun, z, H, fCV = CV, pMut = 0.1)

```

By default, a PMX Crossover and a random swap mutation will be used if the genotype is a Vector{Int}

and a two-point crossover and random flips will be used for Vector{Bool} / BitArrays.

Other crossover / mutations operators can be passed with the keywords `fmut` and `fcross`

```
nsga(popsize, nbGen, init, z , fmut = ..., fcross = ...)
```

Starting solutions are defined with the keyword `seed`
e.g : 
```
x1, x2, x3 = greedy(...)
nsga(popsize, nbGen, init, z, seed = [x1, x2, x3])
```

A plot function can be passed with the keyword `fplot`
e.g. with PyPlot:
```
function plot_pop(P)
    clf()
    p = plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    sleep(0.2)
end

nsga(popsize, nbGen, init, z, fplot = plot_pop)
```

## RealCoding

`RealCoding(eps, lb, ub)` and `decode(x, d::RealCoding)` can be used to easily represent Real values with eps decimal places.

`encode(x, d)` can be used to provide starting solutions

```
d = RealCoding(6, [-3, -3], [6, 6]) #Codes two reals with 6 decimal places with lower bound -3 and upper bound 6
z1(x1, x2) = -(3(1-x1)^2 * exp(-x1^2 - (x2+1)^2) - ...)
z2(x1, x2) = -(3(1+x2)^2 * exp(-x2^2 - (1-x1)^2) - ...)
z(x) = begin 
    x1, x2 = decode(x, d)
    z1(x1, x2), z2(x1, x2)
end
x1, x2 = (-3., 6), (-$\pi$, 0)
seed = encode.([x1, x2], d)
nsga(100, 50, ()->bitrand(d.nbbitstotal), z, seed=seed)
```


## vOptGeneric

This package supports models declared with JuMP / vOptGeneric, with the restriction that all variables must be bounded, and that they're either continuous or binary (integer might get added later)

See examples/Mavrota_MILP.jl
