function DE_evolve(D, fidelity_arr, μ0::Float64, ξ0::Float64, P::Int64)

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
    M = zeros(P, knobs, N)
    for i in 1:P
        r = sample(1:P, 3, replace = false)
        M[i,:,:] = D[r[1],:,:] + μ*(D[r[2],:,:] - D[r[3],:,:])
    end
    println("mutation complete")
    #Crossover
    C = zeros(P, knobs, N)
    for i in 1:P
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
    println("crossover complete")
    #Selection
    #TODO Paralellise this rate determining step !!
    for i in 1:P
        U1 = Integ_H(C[i,:,:])
        U2 = Integ_H(D[i,:,:])
        f1 = fidelity(U1,Utarget)
        f2 = fidelity(U2,Utarget)
        if f1 > f2
            fidelity_arr[i] = f1
            D[i,:,:] = C[i,:,:]
        else
            fidelity_arr[i] = f2
        end
    end
    println("selection complete")
    return μ, ξ
end

function DE_iter()
    D = rand(DE_population, knobs, N)
    fidelity_arr = zeros(DE_population)
    fidelity_ts = zeros(generations)
    μ0 = 0.9
    ξ0 = 0.5

    for i in 1:generations
        println(μ0)
        μ0, ξ0 = DE_evolve(D, fidelity_arr, μ0, ξ0, DE_population)
        fidelity_ts[i] = maximum(fidelity_arr)
    end

    @save "fidelity_ts.jld2" fidelity_ts
    pygui(true)
    plot(fidelity_ts)
end
