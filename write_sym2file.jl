open("A2.txt", "w") do f
    for i in 1:20
       for j in 1:20
           write(f, string(A2_proj[i,j]))
           write(f, " ")
           end
       write(f, ";\n")
       end
    end
