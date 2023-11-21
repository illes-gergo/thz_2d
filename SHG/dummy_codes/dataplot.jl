using HDF5, PlotlyJS

e0 = 8.8541878128e-12
c0 = 3e8

include("fuggvenyek.jl")
include("plotfuncs.jl")

calcdir = "/home/illesg/cst/2d-calculations"
calcdir1D = "/home/illesg/cst/calculations/2.0-ps"

FID = h5open(calcdir * "/100gw2psfull8mm.hdf5", "r")
#FID1D = h5open(calcdir1D * "DB-2ps-100GW", "r")

keys(FID)

lout_general = Layout(xaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), yaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), font=attr(size=14, color="Black"), plot_bgcolor=:white)



En0 = sum(read(FID["/En"])[1, :])
EnTHz = sum(read(FID["/EnTHz"]), dims=2)

effic = EnTHz ./ En0

effic = dropdims(effic, dims=2) * 100

lout_effic = Layout(yaxis=attr(range=[0, 3.5], title="Efficiency (%)"), xaxis=attr(title="Crystal length (mm)"))

#p = Plot(scatter(x=(0:400)*20e-3,y=effic * 4 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) .^ 2),merge(lout_general,lout_effic))

crylength = "175"

ETHz = real(read(FID[crylength*"/ATHz_xt"])) * 2 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) / 1e5

x = read(FID["/x"])
t = read(FID["/t"])

_, I = findmax(abs.(read(FID[crylength*"/ATHz_xt"])))

lout_thz = Layout(xaxis=attr(range=[-3, 3], title="Time (ps)"), yaxis=attr(range=[-7, 10], title="Transversal pos. (mm)"))

lout_thz_spectr = Layout(xaxis=attr(range=[0, 1], title="Frequency (THz)"), yaxis=attr(range=[-7,10], title="Transversal pos. (mm)"))

lout_thz1D = Layout(xaxis=attr(range=[-5, 5], title="Time (ps)"), yaxis=attr(title="Electric field (kV/cm)"), legend=attr(x=0.8, y=1))

lout_thz1D_spectr = Layout(xaxis=attr(range=[0, 5], title="Frequency (THz)"), yaxis=attr(title="Spectral Amplitude (arb. u.)",range=[0,14e-6]), legend=attr(x=0.8, y=1))

lout_imax=Layout(xaxis=attr(range=[0,2],title="Kristályhossz (mm)"),yaxis=attr(title=L"Intenzitás (GW/cm2$)",range=[70,105]))

#= p = Plot(heatmap(x=t * 1e12, y=x * 1e3, z=circshift(ETHz, (512 - I[1], 0))', colorscale="Jet", colorbar_title = "E (kV/cm)"), merge(lout_general, lout_thz)) =#

E_shift = circshift(ETHz, (512 - I[1], 1024 - I[2]))

#y1 = 

#= data1 = (scatter(x=t*1e12, y=E_shift[:,1024],name="x=0 mm"))
data2 = (scatter(x=t*1e12, y=E_shift[:,1024-82], name="x=-2 mm"))
data3 = (scatter(x=t*1e12, y=E_shift[:,1024+82],name="x=-2 mm"))

p = Plot([data1;data2;data3],merge(lout_general,lout_thz1D)) =#

ATHz = read(FID[crylength*"/ATHz_xo"])

freq = read(FID["/omega"]) / 2 / pi

#= p = Plot(heatmap(x=(freq.-freq[1]) / 1e12, y=x * 1e3, z=abs.(ATHz)', colorscale="Jet", colorbar_title="Spect. Amp. (arb. u.)",colorbar_titleside="right"), merge(lout_general, lout_thz_spectr)) =#

#= ATHz_shift = abs.(circshift(ATHz,(0,1024-I[2])))

data1 = (scatter(x=(freq.-freq[1]) / 1e12, y=ATHz_shift[:,1024],name="x=0 mm"))
data2 = (scatter(x=(freq.-freq[1]) / 1e12, y=ATHz_shift[:,1024-82], name="x=-2 mm"))
data3 = (scatter(x=(freq.-freq[1]) / 1e12, y=ATHz_shift[:,1024+82],name="x=-2 mm"))

p = Plot([data1;data2;data3],merge(lout_general,lout_thz1D_spectr)) =#


z = 1:401

IopM = I_opMax(z)

data = scatter(x=z*20e-3,y=IopM)
p = plot(data,merge(lout_general,lout_imax))

#display(p)

close(FID)
#close(FID1D)