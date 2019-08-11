using SymPy
using PyCall

global const O3 = projector_gen()

function commutator(H1,H2)
    #Evaluate the commutator of given Hamiltonian at times t1 and t2
    Hcom = H1*H2 - H2*H1
    return Hcom
end

function h_symb(var_ind)
    η = 0.2
    g = 0.03
    e1, e2, e3 = symbols("e1$(var_ind) e2$(var_ind) e3$(var_ind)")
    I4 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    H1 = [0 0 0 0; 0 e1 0 0; 0 0 2 * e1 - η 0; 0 0 0 3 * e1 - 3 * η;]
    H2 = [0 0 0 0; 0 e2 0 0; 0 0 2 * e2 - η 0; 0 0 0 3 * e2 - 3 * η;]
    H3 = [0 0 0 0; 0 e3 0 0; 0 0 2 * e3 - η 0; 0 0 0 3 * e3 - 3 * η;]
    X = [0 1 0 0; 1 0 sqrt(2) 0; 0 sqrt(2) 0 sqrt(3); 0 0 sqrt(3) 0]
    Y = [0 -1 0 0; 1 0 -sqrt(2) 0; 0 sqrt(2) 0 -sqrt(3); 0 0 sqrt(3) 0] #should be multiplied by an i

    I8 = kron(I4, I4)

    H = kron(H1, I8) + kron(kron(I4, H2), I4) + kron(I8, H3) + g / 2 * ( kron(kron(X, X), I4) + kron(I4, kron(X, X)) - kron(kron(Y, Y), Y) - kron(I4, kron(Y, Y)) )
    return e1,e2,e3,H
end

e11,e21,e31,H1 = h_symb(1)
e12,e22,e32,H2 = h_symb(2)
e13,e23,e33,H3 = h_symb(3)

A1 = H1
A1_O3 = O3*A1*O3'
A1_proj = A1_O3[1:20,1:20]

for i in 1:20
    for j in 1:20
        A1_proj[i,j] = integrate(A1_proj[i,j])
    end
end

A2 = H1*H2 - H2*H1
A2_O3 = O3*A2*O3'
A2_proj = A2_O3[1:20,1:20]

A3 = commutator(H3, commutator(H2,H1)) + commutator(commutator(H3,H2), H1)
A3_O3 = O3*A3*O3'
A3_proj = A3_O3[1:20,1:20]

#=
e12,e22,e32,H2 = h_symb(2)
#e13,e23,e33,H3 = h_symb(3)

A1 = H1
println("clear")
A2 = commutator(H1,H2)
println("clear")
#A3 = commutator(H3, commutator(H2,H1)) + commutator(commutator(H3,H2), H1)
#println("clear")

f_A1 = lambdify(A1,(e11,e21,e31))
f_A1_sparse = lambdify(H1_sparse,(e11,e21,e31))
println("clear")
f_A2 = lambdify(A2,(e11,e21,e31,e12,e22,e32))
println("clear")
#f_A3 = lambdify(A3,(e11,e21,e31,e12,e22,e32,e13,e23,e33))
#println("clear")
=#
