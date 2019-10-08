using PyPlot, JLD2, SharedArrays, Distributed
#@load "/home/vamsi/Github/Quantum-Gate-Implementation-master/data2019-09-19T19:50:38.173.jld2"
f=jldopen("/home/vamsi/Github/Quantum-Gate-Implementation-master/data/99%/data2019-09-17T16:51:17.605.jld2","r")
fidelity_ts = f["fidelity_ts"]
D = f["D"]
pygui(true)
#plot(D[1,:])
figure()
plot(fidelity_ts)

ylabel("Fidelity")
xlabel("Iterations")

figure()
for i in 1:size(D)[1]
    plot(D[i,1,:])
    plot(D[i,2,:])
    plot(D[i,3,:])
end
xlabel("Time(ns)")
ylabel("Control pulse (Ghz)")
