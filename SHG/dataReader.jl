using Plots, JLD2, FFTW

plotlyjs()

db1 = jldopen("23-05-09 10-29-13.jld2")
db2 = jldopen("23-05-09 10-29-12noSHG.jld2")
AxtSHG = db1["/Axt"]
Axt = db2["/Axt"]

p1 = heatmap(abs.(Axt))
p2 = heatmap(abs.(AxtSHG))

M1 = maximum(abs.(Axt))
M2 = maximum(abs.(AxtSHG))

plot(p1,p2)