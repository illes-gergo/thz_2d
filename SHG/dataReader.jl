using Plots, HDF5, FFTW, DelimitedFiles

include("fuggvenyek.jl")

plotlyjs()

FID = h5open("23-05-22 15-00-50.hdf5", "r")

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)

EnergySH = Array{Float64}(undef, 1)
EnergyTHz = Array{Float64}(undef, 1)

maxEntry = read(FID["/maxEntry"])

for i in 1:maxEntry
    push!(EnergySH, sum(abs.(read(FID["/"*string(i)*"/ASH"])) .^ 2))
    push!(EnergyTHz, sum(abs.(read(FID["/"*string(i)*"/ATHz_xo"])) .^ 2))
end


z = range(start=0, stop=8e-3, length=maxEntry + 1)

EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(0.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100

gamma = deg2rad(read(FID["/gamma"]))

display(plot(z, EfficTHz))
FILE = readdlm("efficSH.txt")
plot(z, EfficSH)
display(plot!(FILE[:, 1]/cos(gamma), FILE[:, 2]*50))




close(FID)