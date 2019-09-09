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
    if rand() > S
        M = zeros(DE_population, knobs, N)
        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,:] = D[r[1],:,:] + μ*(D[r[2],:,:] - D[r[3],:,:])
        end
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
    else
        M = zeros(DE_population, knobs, N)
        ctrpar_index = sample(1:N,1)
        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,ctrpar_index] = D[r[1],:,ctrpar_index] + μ*(D[r[2],:,ctrpar_index] - D[r[3],:,ctrpar_index])
        end
        #Crossover
        C = deepcopy(D)
        for i in 1:DE_population
            for k in 1:knobs
                if rand() < ξ
                    C[i,k,ctrpar_index] = M[i,k,ctrpar_index]
                else
                    C[i,k,ctrpar_index] = D[i,k,ctrpar_index]
                end
            end
        end
    end

    #Selection
    #TODO Paralellise this rate determining step !!
    @distributed for i in 1:DE_population
        U1 = Integ_H(C[i,:,:])
        f1 = fidelity(U1,Utarget)
        if f1 > fidelity_arr[i]
            fidelity_arr[i] = f1
            D[i,:,:] = C[i,:,:]
        end
    end

    return μ, ξ
end

function DE_iter()
    D = SharedArray{Float64,3}(DE_population, knobs, N)
    for i in 1:DE_population,j in 1:knobs,k in 1:N
        D[i,j,k] = rand()
    end
    fidelity_arr = SharedArray{Float64}(DE_population)
    fidelity_ts = SharedArray{Float64}(generations)
    μ0 = 0.9
    ξ0 = 0.5

    @distributed for i in 1:DE_population
        U = Integ_H(D[i,:,:])
        fidelity_arr[i] = fidelity(U,Utarget)
    end

    for i in 1:generations
        println(i)
        println(μ0)
        μ0, ξ0 = DE_evolve(D, fidelity_arr, μ0, ξ0)
        fidelity_ts[i] = maximum(fidelity_arr)
    end

    @save "fidelity_ts"*string(Dates.now())*".jld2" fidelity_ts D
    #pygui(true)
    #plot(fidelity_ts)
end
