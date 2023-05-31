using Plots, HDF5, FFTW, DelimitedFiles

include("fuggvenyek.jl")

plotlyjs()

#FID = h5open("23-05-24 11-03-45-no4511G.hdf5", "r")
FID = h5open("23-05-26 13-19-41.hdf5", "r")

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)

EnergySH = Array{Float64}(undef, 1)
EnergyTHz = Array{Float64}(undef, 1)

maxEntry = read(FID["/maxEntry"])

for i in 2:maxEntry
    push!(EnergySH, sum(abs.(read(FID["/"*string(i)*"/ASH"])) .^ 2))
    push!(EnergyTHz, sum(abs.(read(FID["/"*string(i)*"/ATHz_xo"])) .^ 2))
end


z = range(start=0, stop=8e-3, length=maxEntry)

EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(0.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100

gamma = deg2rad(read(FID["/gamma"]))

display(plot(z, EfficTHz))
FILE = readdlm("SHG/efficSH.txt")
plot(z, EfficSH)
display(plot!(FILE[:, 1], FILE[:, 2]*100))

mult = 2.8 

Aop = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/Aop"])
Eop = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/Eop"])

ATHz = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/ATHz_xo"])
ETHz = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/ATHz_xt"])

p1 = (heatmap(abs.(Aop)))
p2 = (heatmap(abs.(Eop)))

p3 = (heatmap(abs.(ATHz)))
p4 = (heatmap(real.(ETHz)))

display(plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900)))

p11 = (plot(abs.(Aop)[:,125]))
p21 = (plot(abs.(Eop)[:,125]))

p31 = (plot(abs.(ATHz)[:,125]))
p41 = (plot(real.(ETHz)[:,125]))

display(plot(p11, p21, p31, p41, layout=(2, 2), size=(1200, 900)))

close(FID)