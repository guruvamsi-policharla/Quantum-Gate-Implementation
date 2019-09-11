using Distributed
addprocs(4)
println(nprocs())
@everywhere using LinearAlgebra, SpecialFunctions, StatsBase, JLD2, Dates, SharedArrays, Statistics

#@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/data/fidelity_ts2019-09-11T22:35:54.6060.14_80_1500.jld2"
@everywhere include("hamiltonian.jl")
@everywhere include("EvalUT.jl")
@everywhere include("DE.jl")

@everywhere global const N = 27
@everywhere global const knobs = 3
@everywhere global const μl = 0.1
@everywhere global const μu = 0.1
@everywhere global const κ1 = 0.1
@everywhere global const κ2 = 0.9
@everywhere global const run_time = 27 #ns
@everywhere global const O3 = projector_gen()
@everywhere global const comp_dim = 20
@everywhere global const DE_population = 20
@everywhere global const generations = 5
@everywhere global const Utarget = [
1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0;
0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
]
@everywhere global const S = 0.0


fidelity_arr = SharedArray{Float64,1}(DE_population)
U = rand(DE_population,20,20)# + im*rand(DE_population,20,20)

@sync @distributed for i in 1:DE_population
    fidelity_arr[i] = fidelity(U[i,:,:],Utarget)
    if i == DE_population
        display(U[i,1:8,1:8])
        println("fidelity",fidelity(U[i,:,:],Utarget))
        display(Utarget)
    end
end

fidelity_arr_temp = zeros(DE_population)
for i in 1:DE_population
    fidelity_arr_temp[i] = fidelity(U[i,:,:],Utarget)
    if i == DE_population
        display(U[i,1:8,1:8])
        println("fidelity",fidelity(U[i,:,:],Utarget))
        display(Utarget)
    end
end

#println(fidelity_arr-fidelity_arr_temp)
rmprocs(workers())
