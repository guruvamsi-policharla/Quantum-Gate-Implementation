include("hamiltonian.jl")
include("EvalUT.jl")

using CSV, DataFrames, PyCall, PhysicalConstants, LinearAlgebra, SpecialFunctions, StatsBase

global const N = 27
global const knobs = 3
global const μl = 0.1
global const μu = 0.1
global const κ1 = 0.1
global const κ2 = 0.9
global const run_time = 27 #ns
global const ħ = 1.0545718176461565e-16 # fix units
global const O3 = projector_gen()
global const comp_dim = 20
