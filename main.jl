using Distributed
addprocs(4)
println(nprocs())
@everywhere using LinearAlgebra, SpecialFunctions, StatsBase, JLD2, Dates, SharedArrays, Statistics, Optim

@everywhere include("hamiltonian.jl")
@everywhere include("EvalUT.jl")
@everywhere include("DE.jl")

@everywhere global const N = 30
@everywhere global const knobs = 3
@everywhere global const μl = 0.1
@everywhere global const μu = 0.1
@everywhere global const κ1 = 0.1
@everywhere global const κ2 = 0.9
@everywhere global const run_time = 30 #ns
@everywhere global const O3 = projector_gen()
@everywhere global const comp_dim = 20
@everywhere global const DE_population = 10
@everywhere global const generations = 5
@everywhere global const Utarget = [
1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
0.0  0.0  0.0  0.0  0.0  0.0  0.0  -1.0;
]
@everywhere global const S = 0.14

#=
Fixing the units as follows
ϵ is in Ghz
T is in ns
Hence, do not need to change h if we ignore the 10^9 factor in both time and frequency
=#

@time DE_iter()

rmprocs(workers())
println("Successfully removed workers")
