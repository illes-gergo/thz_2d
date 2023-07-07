using Plots, HDF5, FFTW, DelimitedFiles, LazyGrids

include("fuggvenyek.jl")

gr()

#FID = h5open("23-05-24 11-03-45-no4511G.hdf5", "r")
FID = h5open("SHG/próbaszámolás.hdf5", "r")

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)

EnergySH = Array{Float64}(undef, 1)
EnergyTHz = Array{Float64}(undef, 1)

maxEntry = read(FID["/maxEntry"])

#= for i in 2:maxEntry
    push!(EnergySH, sum(abs.(read(FID["/"*string(i)*"/ASH"])) .^ 2))
    push!(EnergyTHz, sum(abs.(read(FID["/"*string(i)*"/ATHz_xo"])) .^ 2))
end
 =#

z = range(start=0, stop=8e-3, length=maxEntry)

#= EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(0.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
 =#
gamma = deg2rad(read(FID["/gamma"]))

#= display(plot(z, EfficTHz))
FILE = readdlm("SHG/efficSH.txt")
plot(z, EfficSH)
display(plot!(FILE[:, 1], FILE[:, 2]*100))
 =#
mult = 3

t = read(FID["/t"])
x = read(FID["/x"])

cx, ct = ndgrid(x,t)

Aop = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/Aop"])
Eop = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/Eop"])

ATHz = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/ATHz_xo"])
ETHz = read(FID["/"*string(floor(Int, (maxEntry - 1) / 8 * mult))*"/ATHz_xt"])
#=
p1 = (heatmap(abs.(Aop)))
p2 = (heatmap(abs.(Eop)))

p3 = (heatmap(abs.(ATHz)))
p4 = (heatmap(real.(ETHz)))

display(plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900)))

p11 = (plot(abs.(Aop)[:,1020]))
p21 = (plot(abs.(Eop)[:,1020]))

p31 = (plot(abs.(ATHz)[:,1020]))
p41 = (plot(real.(ETHz)[:,1020]))

display(plot(p11, p21, p31, p41, layout=(2, 2), size=(1200, 900)))
 =#
c0 = 3e8
e0 = 8.854187817e-12
close(FID)
shift = 82
ETHz2 = real.(circshift(ETHz,(shift,0)))'
Eop2 = abs.(circshift(Eop,(shift,0)))'
pltTHz = contour(t[300:500]*1e12,x[800:1201]*1e3,ETHz2[800:1201,300:500]/1e5,xlabel="Time (ps)", ylabel="x (mm)", colorbar_title="\nElectric field (kV/cm)", 
framestyle=:box, colormap=:jet, tickfont=14, guidefont=16, size=(400,300).*2,fill=true,levels=100,lw=0,clims=(-maximum(abs.(ETHz2)),maximum(abs.(ETHz2)))./1e5, 
right_margin = 10Plots.mm, colorbar_titlefontsize=16)
display(pltTHz)

pltOp = contour(t[300:500]*1e12,x[800:1201]*1e3,Eop2[800:1201,300:500].^2/1e13/2*c0*e0*neo(10.6e-6,300,4),xlabel="Time (ps)", ylabel="x (mm)", colorbar_title="\nIntensity (GW/cm2)",
 framestyle=:box, colormap=:jet, tickfont=14, guidefont=16, size=(400,300).*2,fill=true,levels=100,lw=0,right_margin = 10Plots.mm,
 clims=(0,maximum(abs.(Eop2))).^2 ./1e13 ./2 .*c0 .*e0 .*neo(10.6e-6,300,4),colorbar_titlefontsize=16)
 display(pltOp)
#writedlm("tyeragec-cont.txt",[ct[:];;cx[:];;(real(ETHz)')[:]])