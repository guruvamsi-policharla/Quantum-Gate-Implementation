function DE_evolve(μ0, ξ0, P::Int64)
    D = zeros(P, knobs, N)

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

    #Crossover
    C = zeros(P, knobs, N)
    for i in 1:P
        for j in 1:N
            for k in 1:knobs
                if rand() < ξ
                    C[i,knobs,k] = M[i,knobs,j]
                else
                    C[i,knobs,k] = D[i,knobs,j]
                end
            end
        end
    end

    #Selection
    for i in 1:P
        if fidelity(C[i,:,:],Utarget) > fidelity(D[i,:,:],Utarget)
            D[i,:,:] = C[i,:,:]
        #else do nothing
        end
    end

end
