#Working under the assumption that the Hamiltonian is available at all instance of time.

function commutator(H1,H2)
    #Evaluate the commutator of given Hamiltonian at times t1 and t2
    Hcom = H1*H2 - H2*H1
    return Hcom
end

function magnus_exp(H,T)

    N2 = N^2
    dt = T/N

    A1 = zeros(20,20)
    A2 = zeros(20,20)
    A3 = zeros(20,20)

    for i in 0:N-1
        A1 = A1 + H(i*dt)
    end

    A1 = A1*dt

    for i in 0:N-1
        for j in 0:i
            A2 = A2 + H(j*dt)*H(i*dt) - H(i*dt)*H(j*dt)
        end
    end
    A2 = A2*dt^2
    #=
    for i in 0:N-1
        for j in 0:i
            for k in 0:j
                A3 = A3 + commutator( H(k*dt), commutator(H(j*dt), H(i*dt)) ) + commutator( commutator(H(k*dt),H(j*dt)), H(i*dt) )
            end
        end
        if mod(i,10) == 0
            println(i)
        end
    end

    A3 = A3*dt^3
    =#
    return A1*-im/ħ + A2/ħ^2/2# + A3*im/ħ^3/6
end
