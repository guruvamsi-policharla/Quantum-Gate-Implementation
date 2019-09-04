function DE_evolve(D, fidelity_arr, μ0::Float64, ξ0::Float64)

    if rand() < κ1
        μ = μl + rand()*μu
    else
        μ = μ0
    end

    if rand() < κ2
        ξ = rand()
    else
        ξ = ξ0
    end

    #Mutation
    M = zeros(DE_population, knobs, N)
    for i in 1:DE_population
        r = sample(1:DE_population, 3, replace = false)
        M[i,:,:] = D[r[1],:,:] + μ*(D[r[2],:,:] - D[r[3],:,:])
    end
    #println("mutation complete")
    #Crossover
    C = zeros(DE_population, knobs, N)
    for i in 1:DE_population
        for j in 1:N
            for k in 1:knobs
                if rand() < ξ
                    C[i,k,j] = M[i,k,j]
                else
                    C[i,k,j] = D[i,k,j]
                end
            end
        end
    end
    #println("crossover complete")
    #Selection
    #TODO Paralellise this rate determining step !!
    #TODO Optimise by remembering the fidelity
    for i in 1:DE_population
        U1 = Integ_H(C[i,:,:])
        f1 = fidelity(U1,Utarget)
        if f1 > fidelity_arr[i]
            fidelity_arr[i] = f1
            D[i,:,:] = C[i,:,:]
        end
    end
    #println("selection complete")
    return μ, ξ
end

function DE_iter()
    D = rand(DE_population, knobs, N)
    fidelity_arr = zeros(DE_population)
    fidelity_ts = zeros(generations)
    μ0 = 0.9
    ξ0 = 0.5

    for i in 1:DE_population
        U = Integ_H(D[i,:,:])
        fidelity_arr[i] = fidelity(U,Utarget)
    end

    for i in 1:generations
        println(i)
        println(μ0)
        μ0, ξ0 = DE_evolve(D, fidelity_arr, μ0, ξ0)
        fidelity_ts[i] = maximum(fidelity_arr)
    end

    @save "fidelity_ts.jld2" fidelity_ts D
    pygui(true)
    plot(fidelity_ts)
end
