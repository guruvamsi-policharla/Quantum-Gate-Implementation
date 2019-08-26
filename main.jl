include("hamiltonian.jl")
include("EvalUT.jl")

using CSV, DataFrames, PyCall, PhysicalConstants, LinearAlgebra, SpecialFunctions

global const N = 27
global ϵ1_nodes = rand(N)
global ϵ2_nodes = rand(N)
global ϵ3_nodes = rand(N)
global const run_time = 27 #ns
global const ħ = 1.0545718176461565e-16 # fix units
global const O3 = projector_gen()
global const comp_dim = 20
