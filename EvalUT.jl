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
    return Utot.*sqrt(comp_dim)./sqrt(tr(Utot'*Utot))
end

function infidelity(U1,U2,x)
    p1 = Diagonal([1, exp(-im*x[6]), exp(-im*x[5]), exp(-im*(x[5] + x[6])), exp(-im*x[4]), exp(-im*(x[4]+x[6])), exp(-im*(x[4]+x[5])), exp(-im*(x[4]+x[5]+x[6]))])
    p2 = Diagonal([1, exp(-im*x[3]), exp(-im*x[2]), exp(-im*(x[2] + x[3])), exp(-im*x[1]), exp(-im*(x[1]+x[3])), exp(-im*(x[1]+x[2])), exp(-im*(x[1]+x[2]+x[3]))])

    _u1 = U1[1:8,1:8]
    _u2 = U2[1:8,1:8]
    return 1-1/8 * abs(tr(p1*_u1*p2 * _u2'))
end

function phasecomp_infidel(U1,U2)
    res = optimize(x -> infidelity(U1,U2,x),[0.1 0.2 0.3 0.1 0.2 0.3])
    return res.minimum
end
