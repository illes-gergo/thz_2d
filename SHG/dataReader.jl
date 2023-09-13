using PlotlyJS, HDF5, FFTW, DelimitedFiles, LazyGrids

include("fuggvenyek.jl")


#FID = h5open("23-05-24 11-03-45-no4511G.hdf5", "r")
FID = h5open("SHG/próbaszámolás.hdf5", "r")

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)


EnergySH = Array{Float64}(undef, 1)
EnergyTHz = Array{Float64}(undef, 1)
EnergyPump = Array{Float64}(undef, 1)
EnergyPump[1] = 0
EnergyTHz[1] = 0
EnergySH[1] = 0

maxEntry = read(FID["/maxEntry"])

for i in 2:maxEntry
    push!(EnergySH, sum(abs.(read(FID["/"*string(i)*"/ASH"])) .^ 2))
    push!(EnergyTHz, sum(abs.(read(FID["/"*string(i)*"/ATHz_xo"])) .^ 2))
    push!(EnergyPump, sum(abs.(read(FID["/"*string(i)*"/Aop"])) .^ 2))
end


z = range(start=0, stop=2e-3, length=maxEntry)

EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(0.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100

gamma = deg2rad(read(FID["/gamma"]))

efficSH1D = readdlm("efficSH1D.txt")

a = plot(scatter(x=z,y=EfficTHz))
b = plot([scatter(x=z,y=EfficSH,name="2D-s SHG hatásfok");scatter(x=efficSH1D[:,1],y=efficSH1D[:,2]*100,name="1D-s SHG hatásfok")])
c = plot(scatter(x=z,y=EnergyPump))
display(b)
display(c)
mult = 3

t = read(FID["/t"])
x = read(FID["/x"])

cx, ct = ndgrid(x, t)

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
ETHz2 = real.(circshift(ETHz, (shift, 0)))'
Eop2 = abs.(circshift(Eop, (shift, 0)))'



#writedlm("tyeragec-cont.txt",[ct[:];;cx[:];;(real(ETHz)')[:]])