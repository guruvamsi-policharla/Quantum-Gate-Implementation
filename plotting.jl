using PyPlot, JLD2

@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/fidelity_ts2019-09-05T20:58:58.647.jld2"
pygui(true)
#plot(D[1,:])
#figure()
plot(fidelity_ts)

ylabel("Fidelity")
xlabel("Iterations")
#=
figure()
for i in 1:5
    plot(D[i,2,:])
end
=#
