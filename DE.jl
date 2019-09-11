function DE_evolve(D::SharedArray, fidelity_arr::SharedArray, μ0::Float64, ξ0::Float64)
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

    M = SharedArray{Float64,3}(DE_population, knobs, N)
    C = SharedArray{Float64,3}(DE_population, knobs, N)

    #Mutation
    if rand() > S
        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,:] = D[r[1],:,:] + μ*(D[r[2],:,:] - D[r[3],:,:])
        end
        #Crossover
        for i in 1:DE_population
            for j in 1:N
                if rand() < ξ
                    for k in 1:knobs
                        C[i,k,j] = M[i,k,j]
                    end
                else
                    for k in 1:knobs
                        C[i,k,j] = D[i,k,j]
                    end
                end
            end
        end

        #Selection
        #TODO Paralellise this rate determining step !!
        @sync @distributed for i in 1:DE_population
            U1 = Integ_H(C[i,:,:])
            f1 = fidelity(U1,Utarget)

            if f1 > fidelity_arr[i]
                fidelity_arr[i] = f1
                D[i,:,:] = C[i,:,:]
            end
        end
    else
        ctrpar_index = sample(1:N,1)

        for i in 1:DE_population,j in 1:knobs,k in 1:N
            M[i,j,k] = D[i,j,k]
            C[i,j,k] = D[i,j,k]
        end

        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,ctrpar_index] = D[r[1],:,ctrpar_index] + μ*(D[r[2],:,ctrpar_index] - D[r[3],:,ctrpar_index])
        end
        #Crossover
        for i in 1:DE_population
                if rand() < ξ
                    for k in 1:knobs
                        C[i,k,ctrpar_index] = M[i,k,ctrpar_index]
                    end
                else
                    for k in 1:knobs
                        C[i,k,ctrpar_index] = D[i,k,ctrpar_index]
                    end
                end
        end

        #Selection
        #TODO Paralellise this rate determining step !!
        @sync @distributed for i in 1:DE_population
            U1 = Integ_H(C[i,:,:])
            f1 = fidelity(U1,Utarget)

            if f1 > fidelity_arr[i]
                fidelity_arr[i] = f1
                D[i,:,:] = C[i,:,:]
            end
        end
    end

    return μ, ξ, D
end

function DE_iter()
    D = SharedArray{Float64,3}(DE_population, knobs, N)
    for i in 1:DE_population,j in 1:knobs,k in 1:N
        D[i,j,k] = rand() - 0.5
    end
    fidelity_arr = SharedArray{Float64,1}(DE_population)
    fidelity_ts = SharedArray{Float64,1}(generations)
    μ0 = 0.9
    ξ0 = 0.5


    @sync @distributed for i in 1:DE_population
        U = Integ_H(D[i,:,:])
        fidelity_arr[i] = fidelity(U,Utarget)
    end

    println("Setup Complete. Evolution Starting.")

    for i in 1:generations
        println(i)
        μ0, ξ0, D = DE_evolve(D, fidelity_arr, μ0, ξ0)
        fidelity_ts[i] = maximum(fidelity_arr)
        println(fidelity_ts[i])
    end
    println("Evolution Complete")

    rmprocs(workers())

    @save "fidelity_ts"*string(Dates.now())*".jld2" fidelity_ts D

end
