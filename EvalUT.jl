function EvalUT(T)
    mesh = 10
    Utot = Matrix{Float64}(I, comp_dim, comp_dim)
    dt = T/N/mesh

    for par_itr in 1:N-1 #iterating through the ctrpars
        #TODO define the interpolation function
        interp_arr_1 = interp(ϵ1_nodes[par_itr], ϵ1_nodes[par_itr+1])
        interp_arr_2 = interp(ϵ2_nodes[par_itr], ϵ2_nodes[par_itr+1])
        interp_arr_3 = interp(ϵ3_nodes[par_itr], ϵ3_nodes[par_itr+1])
        for ϵ_itr in 1:mesh
            h_temp = hamiltonian_o3_eval(interp_arr_1[ϵ_itr], interp_arr_2[ϵ_itr], interp_arr_3[ϵ_itr])
            Utot = Utot * exp(-1*im*dt*h_temp)
        end
    end

    return Utot
end
