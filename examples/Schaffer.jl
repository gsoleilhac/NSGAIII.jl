    # David Schaffer. Multiple Objective Optimization with Vector Evaluated Genetic Algorithms.
    # In Proceedings of the 1st International Conference on Genetic Algorithms, 
    # L. Erlbaum Associates Inc. pp. 93–100, 1985.
    # http://dl.acm.org/citation.cfm?id=645511.657079

    using NSGAIII

    # using PyPlot
    # function plot_pop(P)
    #     clf()
    #     p = plot(map(x -> x.y[1], P), map(x -> x.y[2], P), "bo", markersize=1)
    #     !isinteractive() && show()
    #     sleep(0.2)
    # end

    using UnicodePlots
    function plot_pop(P)
        println()
        display(scatterplot(map(x -> x.y[1], P), map(x -> x.y[2], P)))
        sleep(0.4)
    end

    const d = RealCoding(6, [-10], [10])
    z1(x) = x^2
    z2(x) = (x-2)^2
    z(bits) = begin 
        x = decode(bits, d)[1]
        z1(x), z2(x)
    end

    pop = nsga(100, 20, ()->rand(Bool, d.nbbitstotal), z, 2, fplot = plot_pop)
    pop = nsga(100, 20, ()->rand(Bool, d.nbbitstotal), z, [[1.,0.],[0.5,0.5],[0.,1.]], fplot = plot_pop)

    pop = nsga(100, 20, ()->rand(Bool, d.nbbitstotal), z, 30, fplot = plot_pop)

    [(decode(ind.x, d)[1], ind.y) for ind in pop]
