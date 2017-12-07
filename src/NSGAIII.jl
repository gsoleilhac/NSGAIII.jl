module NSGAIII

using ProgressMeter
using StaticArrays

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("realcoding.jl")

export nsga, RealCoding, decode, encode, DennisDas


function nsga(popSize::Integer, nbGen::Integer, init::Function, z::Function ;
 references=DennisDas(length(z(init())), 3), fCV=(x)->0, pmut=0.05, fmut=default_mutation!, fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)

    @assert popSize % 4 == 0    

    X = typeof(init())
    P = [indiv(init(), z) for _=1:popSize-length(seed)]
    append!(P, indiv.(convert.(X, seed),z))
    fast_non_dominated_sort!(P)
    Q = similar(P)

    @showprogress 0.1 for gen = 1:nbGen

        ind_Q = 1
        for pass = 1:2
            pass == 2 && shuffle!(P)
            for i = 1:4:popSize
                pa = niching_based_selection(P[i], P[i+1], references)
                pb = niching_based_selection(P[i+2], P[i+3], references)
                ca,cb = crossover(pa, pb, fcross)

                rand() < pmut && mutate!(ca, fmut)
                rand() < pmut && mutate!(cb, fmut)

                eval!(ca, z)
                eval!(cb, z)

                Q[ind_Q] = ca
                Q[ind_Q+1] = cb

                ind_Q += 2
            end
        end

        @assert length(Q) == length(P)

        F = fast_non_dominated_sort!(vcat(P, Q))
        i = 1
        empty!(P)
        while length(P) + length(F[i]) <= popSize
            append!(P, F[i])
            i += 1
        end
        
        T = vcat(P, F[i])
        if length(P) != popSize
            normalize_pop!(T)
            associate_references!(T, references)
            n_count = niche_count(P, references)
            while length(P) < popSize
                niching!(P, T, n_count)
            end
        end

        fplot(P)
    end
    P
end

end # module
