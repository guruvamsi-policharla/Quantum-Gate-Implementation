using PyPlot, JLD2, SharedArrays, Distributed
@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/data/data2019-09-12T05:11:09.173.jld2"
#f=jldopen("/home/vamsi/Github/Quantum-Gate-Implementation-master/data/fidelity_ts_5000itr_100pop.jld2","r")
#fidelity_ts = f["fidelity_ts"].s
#D = f["D"].s
pygui(true)
#plot(D[1,:])
figure()
plot(fidelity_ts)

ylabel("Fidelity")
xlabel("Iterations")

figure()
for i in 1:size(D)[1]
    plot(D[i,2,:])
end
