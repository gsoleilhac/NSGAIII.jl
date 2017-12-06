function fast_non_dominated_sort!(pop::Vector{T}) where {T}
    F = T[]
    for p in pop
        empty!(p.dom_list)
        p.dom_count = 0
        for q in pop
            if p ⋖ q
                push!(p.dom_list, q)
            elseif q ⋖ p
                p.dom_count += 1
            end
        end
        if p.dom_count == 0
            p.rank = 1
            push!(F, p)
        end
    end
    res = Vector{T}[]
    i = 2
    while !isempty(F)
        Q = T[]
        for p in F
            for q in p.dom_list
                q.dom_count -= 1
                if q.dom_count == 0
                    q.rank = i
                    push!(Q, q)
                end
            end
        end
        i += 1
        push!(res, F)
        F = Q
    end
    res
end


function niching_based_selection(p1, p2, H)
    
    if p1.CV != p2.CV
        return p1.CV < p2.CV ? p1 : p2
    end

    if p1.ref_point != p2.ref_point
        return rand(Bool) ? p1 : p2
    end

    if p1.rank != p2.rank
        return p1.rank < p2.rank ? p1 : p2
    end

    return distance(p1, H[p1.ref_point]) < distance(p2, H[p2.ref_point]) ? p1 : p2
    
end

distance(a,b) = sqrt(sum((a - b) .^ 2))



function DennisDas(nbobj, p)
    vals = [i for i = 0:1/p:1]
    vcurrent = [[v] for v in vals]
    for i = 2:nbobj
        vnext = Vector{Vector{Float64}}()
        for x in vcurrent
            for v in vals
                if sum(x) + v <= 1.
                    push!(vnext, vcat(x, v))
                end
            end
        end
        vcurrent = vnext
    end
    filter(x -> sum(x) == 1., vcurrent)
end


function normalize_pop!(pop::Vector{T}) where T
    maximums = [i for i in pop[1].y]
    minimums = [i for i in pop[1].y]

    for indiv in pop
        for obj = 1:length(maximums)
            maximums[obj] = max(maximums[obj], indiv.y[obj])
            minimums[obj] = min(minimums[obj], indiv.y[obj])
        end
    end

    for i = 1:length(maximums)
        if maximums[i] == minimums[i]
            minimums[i] -= 1
        end
    end

    for indiv in pop
        indiv.normalized_y = [((indiv.y[i] - minimums[i] / maximums[i] - minimums[i]) for i in 1:length(maximums))...]
    end

end


function associate_references(pop, references)
    for ind in pop

    end
end