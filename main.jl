using CSV, DataFrames, PyCall, PhysicalConstants, LinearAlgebra, SpecialFunctions, StatsBase, PyPlot, JLD2

include("hamiltonian.jl")
include("EvalUT.jl")
include("DE.jl")

global const N = 27
global const knobs = 3
global const μl = 0.1
global const μu = 0.1
global const κ1 = 0.1
global const κ2 = 0.9
global const run_time = 27 #ns
global const O3 = projector_gen()
global const comp_dim = 20
global const DE_population = 20
global const generations = 1000
global const Utarget = Matrix{Float64}(I, 20, 20)

#=
Fixing the units as follows
ϵ is in Ghz
T is in ns
Hence, do not need to change h if we ignore the 10^9 factor in both time and frequency
=#

DE_iter()

#@load "fidelity_ts.jld2"
