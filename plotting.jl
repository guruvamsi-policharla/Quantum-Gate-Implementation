using PyPlot, JLD2

@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/data/fidelity_ts2019-09-02T11:01:07.597.jld2"
pygui(true)
#plot(D[1,:])
plot(fidelity_ts)

ylabel("Fidelity")
xlabel("Iterations")
