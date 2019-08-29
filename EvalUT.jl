function interp(a,b,mesh)
#interpolates between a and b and returns an array with mesh number of equally spaced points in between a and b
    dT = run_time/(N-1)
    time_arr = 0:dT/mesh:dT-1/mesh
    interp_arr = zeros(mesh)

    for i in 1 :mesh
        interp_arr[i] = (a+b)/2 + (b-a)/2 * erf(5/dT * (time_arr[i] - dT/2))
    end

    return interp_arr
end

function Integ_H(ϵ_mat)
    mesh = 10
    Utot = Matrix{Float64}(I, comp_dim, comp_dim)
    dt = run_time/N/mesh

    for par_itr in 1:N-1 #iterating through the ctrpars
        interp_arr_1 = interp(ϵ_mat[1,par_itr], ϵ_mat[1,par_itr+1],mesh)
        interp_arr_2 = interp(ϵ_mat[2,par_itr], ϵ_mat[2,par_itr+1],mesh)
        interp_arr_3 = interp(ϵ_mat[3,par_itr], ϵ_mat[3,par_itr+1],mesh)
        for ϵ_itr in 1:mesh
            h_temp = hamiltonian_o3_eval(interp_arr_1[ϵ_itr], interp_arr_2[ϵ_itr], interp_arr_3[ϵ_itr])
            Utot = exp(-1*im*dt*h_temp) * Utot
        end
    end
    return Utot
end

function fidelity(U1,U2)
    return 1/comp_dim * abs(tr(U1' * U2))
end
