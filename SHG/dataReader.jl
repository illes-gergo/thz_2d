using PlotlyJS, HDF5, FFTW, DelimitedFiles, LazyGrids

include("fuggvenyek.jl")

#FID = h5open("23-05-24 11-03-45-no4511G.hdf5", "r")
FID = h5open("próbaszámolás-koherencia-ss.hdf5", "r")
FID1D = h5open("DB-1ps-60GW", "r")

z1D = dropdims(read(FID1D["/z"])',dims=2)
effic1D = read(FID1D["/effic"])

Energy0 = sum(abs.(collect(FID["1/Aop"])) .^ 2)

maxEntry = read(FID["/maxEntry"])'

EnergySH = zeros(maxEntry)
EnergyTHz = zeros(maxEntry)
EnergySH2 = zeros(maxEntry)
EnergySH3 = zeros(maxEntry)
EnergyPump = zeros(maxEntry)

EnergyPump[1] = sum(abs.(read(FID["/"*string(1)*"/Aop"])) .^ 2)


for i in 1:maxEntry
	EnergySH[i] = sum(abs.(read(FID["/"*string(i)*"/ASH"]))[:,1024] .^ 2)
	EnergyTHz[i] = sum(abs.(read(FID["/"*string(i)*"/ATHz_xo"])) .^ 2)
	EnergyPump[i] = sum(abs.(read(FID["/"*string(i)*"/Aop"])) .^ 2)
	
	EnergySH2[i] = sum(abs.(read(FID["/"*string(i)*"/ASH"]))[:,1084] .^ 2)
	EnergySH3[i] = sum(abs.(read(FID["/"*string(i)*"/ASH"]))[:,954] .^ 2)
end


z_end = read(FID["/z"])[end]

z = range(start=0, stop=z_end, length=maxEntry)

p1 = scatter(x=z,y=EnergySH)
p2 = scatter(x=z,y=EnergySH2)
p3 = scatter(x=z,y=EnergySH3)

EfficSH = EnergySH ./ Energy0 .* neo(10.6e-6 / 2, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100
EfficTHz = EnergyTHz ./ Energy0 .* nTHzo(1.5e12 * 2 * pi, 300, 4) ./ neo(10.6e-6, 300, 4) .* 100

#EnergySH = EnergySH .* neo(10.6e-6 / 2, 300, 4) .* 100
#EnergyTHz = EnergyTHz .* nTHzo(1.5e12 * 2 * pi, 300, 4) .* 100
#EnergyPump = EnergyPump .* neo(10.6e-6, 300, 4) .* 100


gamma = deg2rad(read(FID["/gamma"]))

efficSH1D = dropdims(read(FID1D["/efficSH"])',dims=2)

a = plot(scatter(x=z,y=EfficTHz))
b = plot([scatter(x=z,y=EfficSH,name="2D-s SHG hatásfok");scatter(x=z1D,y=efficSH1D*100,name="1D-s SHG hatásfok")])
c = plot(scatter(x=z,y=EnergyPump))
#display(b)
#display(c)
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
#close(FID)
shift = 82
ETHz2 = real.(circshift(ETHz, (shift, 0)))'
Eop2 = abs.(circshift(Eop, (shift, 0)))';



#writedlm("tyeragec-cont.txt",[ct[:];;cx[:];;(real(ETHz)')[:]])
