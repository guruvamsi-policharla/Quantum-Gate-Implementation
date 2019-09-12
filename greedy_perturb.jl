using Distributed
addprocs(4)
println(nprocs())
@everywhere using LinearAlgebra, SpecialFunctions, StatsBase, JLD2, Dates, SharedArrays, Statistics

@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/data2019-09-11T20:30:17.369.jld2"
@everywhere include("hamiltonian.jl")
@everywhere include("EvalUT.jl")
@everywhere include("DE.jl")
#=
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
=#

fidelity_arr = zeros(DE_population)
for i in 1:DE_population
    U = Integ_H(D[i,:,:])
    fidelity_arr[i] = fidelity(U,Utarget)
end



println(maximum(fidelity_arr))
println(argmax(fidelity_arr))


rmprocs(workers())
