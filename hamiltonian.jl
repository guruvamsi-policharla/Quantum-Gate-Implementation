function ϵ1(t)
    index = convert(Int64,floor(t/run_time*N))
    if index == N+1
        index = N
    end
    return ϵ1_nodes[index+1]
end

function ϵ2(t)
    index = convert(Int64,floor(t/run_time*N))
    if index == N+1
        index = N
    end
    return ϵ2_nodes[index+1]
end

function ϵ3(t)
    index = convert(Int64,floor(t/run_time*N))
    if index == N+1
        index = N
    end
    return ϵ3_nodes[index+1]
end


function hamiltonian_eval(t)
    e1 = ϵ1(t)
    e2 = ϵ2(t)
    e3 = ϵ3(t)

    η = 0.2
    g = 0.03
    I4 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    H1 = [0 0 0 0; 0 e1 0 0; 0 0 2 * e1 - η 0; 0 0 0 3 * e1 - 3 * η;]
    H2 = [0 0 0 0; 0 e2 0 0; 0 0 2 * e2 - η 0; 0 0 0 3 * e2 - 3 * η;]
    H3 = [0 0 0 0; 0 e3 0 0; 0 0 2 * e3 - η 0; 0 0 0 3 * e3 - 3 * η;]
    X = [0 1 0 0; 1 0 sqrt(2) 0; 0 sqrt(2) 0 sqrt(3); 0 0 sqrt(3) 0]
    Y = [0 -1 0 0; 1 0 -sqrt(2) 0; 0 sqrt(2) 0 -sqrt(3); 0 0 sqrt(3) 0] #should be multiplied by an i

    I8 = kron(I4, I4)

    H = kron(H1, I8) + kron(kron(I4, H2), I4) + kron(I8, H3) + g / 2 * ( kron(kron(X, X), I4) + kron(I4, kron(X, X)) - kron(kron(Y, Y), Y) - kron(I4, kron(Y, Y)) )

    return H
end

function hamiltonian_o3_eval(t)
    e1 = ϵ1(t)
    e2 = ϵ2(t)
    e3 = ϵ3(t)

    η = 0.2
    g = 0.03
    I4 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    H1 = [0 0 0 0; 0 e1 0 0; 0 0 2 * e1 - η 0; 0 0 0 3 * e1 - 3 * η;]
    H2 = [0 0 0 0; 0 e2 0 0; 0 0 2 * e2 - η 0; 0 0 0 3 * e2 - 3 * η;]
    H3 = [0 0 0 0; 0 e3 0 0; 0 0 2 * e3 - η 0; 0 0 0 3 * e3 - 3 * η;]
    X = [0 1 0 0; 1 0 sqrt(2) 0; 0 sqrt(2) 0 sqrt(3); 0 0 sqrt(3) 0]
    Y = [0 -1 0 0; 1 0 -sqrt(2) 0; 0 sqrt(2) 0 -sqrt(3); 0 0 sqrt(3) 0] #should be multiplied by an i
    I8 = kron(I4, I4)

    H = kron(H1, I8) + kron(kron(I4, H2), I4) + kron(I8, H3) + g / 2 * ( kron(kron(X, X), I4) + kron(I4, kron(X, X)) - kron(kron(Y, Y), Y) - kron(I4, kron(Y, Y)) )

    H = O3*H*O3'
    H = H[1:20,1:20]

    return H
end

function projector_gen()
    basis_vec = zeros(64)
    itr = 1
    for i in 0:3
        for j in 0:3
            for k in 0:3
                basis_vec[itr] = 100 * i + 10 * j + k
                itr = itr + 1
            end
        end
    end

    P = zeros(64, 64)

    completed = 1
    for i in 0:9
        energy = i
        for j in 1:64
            if(floor(basis_vec[j] / 100) + mod(floor(basis_vec[j] / 10), 10) + mod(basis_vec[j], 10) == energy)
                P[completed,j] = 1;
                completed = completed + 1;
            end
        end
    end

    #P[64,64] = 1
    CSV.write("Permutation_matrix.csv",  DataFrame(P), writeheader = false)
    return P
end
