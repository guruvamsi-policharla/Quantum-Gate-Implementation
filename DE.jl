function DE_evolve(D::SharedArray, fidelity_arr::SharedArray, μ, ξ)
    for i in 1:DE_population
        if rand() < κ1
            μ[i] = μl + rand()*μu
        else
            #Do nothing
        end

        if rand() < κ2
            ξ[i] = rand()
        else
            #Do nothing
        end
    end

    M = SharedArray{Float64,3}(DE_population, knobs, N)
    C = SharedArray{Float64,3}(DE_population, knobs, N)

    #Mutation
    if rand() > S
        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,:] = D[r[1],:,:] + μ[i]*(D[r[2],:,:] - D[r[3],:,:])
        end
        #Crossover
        for i in 1:DE_population
            for j in 1:N
                for k in 1:knobs
                    if rand() < ξ[i]
                        C[i,k,j] = M[i,k,j]
                    else
                        C[i,k,j] = D[i,k,j]
                    end
                end
            end
        end
    else
        ctrpar_index = sample(1:N,1)
        for i in 1:DE_population,k in 1:knobs,j in 1:N
            M[i,k,j] = D[i,k,j]
            C[i,k,j] = D[i,k,j]
        end

        for i in 1:DE_population
            r = sample(1:DE_population, 3, replace = false)
            M[i,:,ctrpar_index] = D[r[1],:,ctrpar_index] + μ[i]*(D[r[2],:,ctrpar_index] - D[r[3],:,ctrpar_index])
        end
        #Crossover
        for i in 1:DE_population
            for k in 1:knobs
                if rand() < ξ[i]
                    C[i,k,ctrpar_index] = M[i,k,ctrpar_index]
                else
                    C[i,k,ctrpar_index] = D[i,k,ctrpar_index]
                end
            end
        end
    end

    #Selection
    @sync @distributed for i in 1:DE_population
        U1 = Integ_H(C[i,:,:])
        f1 = 1 - phasecomp_infidel(U1,Utarget)
        if f1 > fidelity_arr[i]
            fidelity_arr[i] = f1
            D[i,:,:] = C[i,:,:]
        end
    end

    return μ, ξ, D
end

function DE_iter()
    D = SharedArray{Float64,3}(DE_population, knobs, N)
    fidelity_arr = SharedArray{Float64,1}(DE_population)
    fidelity_ts = SharedArray{Float64,1}(generations)
    μ = 0.9*ones(DE_population)
    ξ = 0.5*ones(DE_population)

    while(maximum(fidelity_arr) < 0.6)
        for i in 1:DE_population,j in 1:knobs,k in 1:N
            D[i,j,k] = rand() - 0.5
        end

        @sync @distributed for i in 1:DE_population
            U = Integ_H(D[i,:,:])
            fidelity_arr[i] = 1 - phasecomp_infidel(U,Utarget)
        end
        println("Trying setup again")
        println(maximum(fidelity_arr))

    end
    println("Setup Complete. Evolution Starting.")

    for i in 1:generations
        println(i)
        μ, ξ, D = DE_evolve(D, fidelity_arr, μ, ξ)
        fidelity_ts[i] = maximum(fidelity_arr)

        if(maximum(var(D,dims=1)) < 10^-5)
            break
        end

        println(fidelity_ts[i])
        println("Variance-",maximum(var(D,dims=1)))
    end
    println("Evolution Complete")

    _D = D.s
    _fidelity_ts = fidelity_ts.s
    @save "data"*string(Dates.now())*".jld2" _fidelity_ts _D N knobs run_time O3 comp_dim DE_population generations Utarget S
    println("Data saved")
end
