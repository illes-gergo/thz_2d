using Plots, HDF5, FFTW

plotlyjs()

FID = h5open("SHG/23-05-16 17-10-35.hdf5", "r")

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)

EnergySH = Array{Float64}(undef, 1)
EnergyTHz = Array{Float64}(undef, 1)
for i in 1:61
    push!(EnergySH, sum(abs.(collect(FID["/"*string(i)*"/ASH"])) .^ 2))
    push!(EnergyTHz, sum(abs.(collect(FID["/"*string(i)*"ATHz_xo"])) .^ 2))
end

EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(0.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100

display(plot(EfficSH))
display(plot(EfficTHz))

close(FID)