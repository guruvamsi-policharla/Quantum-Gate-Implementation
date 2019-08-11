include("timeOrdering.jl")
include("hamiltonian.jl")
using CSV, DataFrames, PyCall, PhysicalConstants, LinearAlgebra

global const N = 100 #ns
global ϵ1_nodes = rand(N)
global ϵ2_nodes = rand(N)
global ϵ3_nodes = rand(N)
global const run_time = 30 #ns
global const ħ = 1.0545718176461565e-16 # per ns and units of frequency in Ghz
global const O3 = projector_gen()

@time A = magnus_exp(hamiltonian_o3_eval,run_time);
