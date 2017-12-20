module NSGAIII

using ProgressMeter
using Requires

include("indivs.jl")
include("functions.jl")
include("crossover.jl")
include("mutation.jl")
include("realcoding.jl")

export nsga, RealCoding, decode, encode, DennisDas

function nsga(popSize::Integer, nbGen::Integer, init::Function, z::Function, H::Int=5, fCV=(x)->0.; args...)
    nbobj = length(z(init()))
    nsga(popSize, nbGen, init, z, DennisDas(nbobj, H) ,fCV ; args...)
end

function nsga(popSize::Integer, nbGen::Integer, init::Function, z::Function, references::Vector{Vector{Float64}} ,fCV ;
  pmut=0.05, fmut=default_mutation!, fcross = default_crossover, seed=typeof(init())[], fplot = (x)->nothing)

    popSize = max(popSize, length(references)) 
    if popSize % 4 != 0
        popSize += 4 - (popSize % 4)
    end

    X = typeof(init())
    P = [indiv(init(), z, fCV) for _=1:popSize-length(seed)]
    append!(P, indiv.(convert.(X, seed),z, fCV))
    fast_non_dominated_sort!(P)
    associate_references!(P, references)
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

                eval!(ca, z, fCV)
                eval!(cb, z, fCV)

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
        

        if length(P) != popSize
            T = vcat(P, F[i])
            normalize_pop!(T)
            associate_references!(T, references)
            n_count = niche_count(P, references)
            while length(P) < popSize
                niching!(P, F[i], n_count)
            end
        else
            associate_references!(P, references)
        end

        fplot(P)
    end
    P
end



@require vOptGeneric begin
    function nsga(popSize, nbGen, H, m ; kwargs...)

        vd = @eval Main getvOptData(m)
        @assert all(isfinite, m.colLower)
        @assert all(isfinite, m.colUpper)

        @assert !(:Int in m.colCat) "Only continuous and binary variables are supported"
        ϵ = map(x -> x==:Cont ? 4 : 0, m.colCat)
        d = RealCoding(ϵ, m.colLower, m.colUpper)

        init = () -> bitrand(d.nbbitstotal)

        function evaluate(obj, x)
            dot(obj.aff.coeffs, x[map(v-> getfield(v, :col), obj.aff.vars)]) + obj.aff.constant
        end

        function z(bits)
            x = decode(bits, d)
            ((evaluate(obj, x) for obj in vd.objs)...)
        end

        function CV(bits)
            x = decode(bits, d)
            res = 0.
            for CSTR in m.linconstr

                if CSTR.lb != -Inf && CSTR.lb != typemin(Float64) 
                    if CSTR.lb == 0
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)])
                        res += max(0, -g)
                    elseif CSTR.lb > 0 
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.lb - 1
                        res += max(0, -g)
                    else
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.lb - 1
                        res += max(0, g)
                    end

                end

                if CSTR.ub != Inf && CSTR.lb != typemax(Float64)
                     if CSTR.ub == 0
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)])
                        res += max(0, g)
                    elseif CSTR.ub > 0 
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.ub - 1
                        res += max(0, g)
                    else
                        g = dot(CSTR.terms.coeffs, x[map(v-> getfield(v, :col), CSTR.terms.vars)]) / CSTR.ub - 1
                        res += max(0, -g)
                    end
                end
            end
            res
        end


        for i = 1:length(vd.objs)
            if vd.objSenses[i] == :Max
                vd.objs[i] = vd.objs[i] * -1
            end
        end

        res = nsga(popSize, nbGen, init, z, H, CV ; kwargs...)

        for i = 1:length(vd.objs)
            if vd.objSenses[i] == :Max
                vd.objs[i] = vd.objs[i] * -1
            end
        end

        # signs = Tuple(s == :Min ? 1 : -1 for s in vd.objSenses)

        # for indiv in res
        #     indiv.y = indiv.y .* signs
        # end

        [(decode(ind.x, d), ind.y, ind.CV) for ind in res]

    end
end

end # module
