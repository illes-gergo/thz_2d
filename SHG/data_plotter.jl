using PlotlyJS, HDF5, FFTW, Interpolations
include("data_plotter_layouts.jl")

calcdir = "/home/illesg/cst/2d-calculations/"

entries = readdir(calcdir)
#println(entries)

fpath = calcdir * "/100gw2psfull8mm.hdf5"

FID = h5open(fpath, "r")

EnSH = dropdims(sum(read(FID["/EnSH"]), dims=2), dims=2)
En0 = sum(read(FID["/En"])[1, :])
z = read(FID["z"])

efficSH = EnSH ./ En0 * 1e2

itp = interpolate(efficSH, BSpline(Cubic(Line(OnGrid()))))

sitp = scale(itp,range(z[1],z[end],401))

trace = scatter(x=z*1e3, y=sitp.(z))

plt = Plot(trace, merge(l_generic, l_efficsh))

display(plt)

step = 14

Eop1 = abs.(read(FID["$(step)/Eop"])).^2
Eop2 = abs.(read(FID["$(step+5)/Eop"])).^2

_, I1 = findmax(Eop1)
lEop1 = circshift(Eop1[:,I1[2]],(512-I1[1]))

_, I2 = findmax(Eop2)
lEop2 = circshift(Eop2[:,I2[2]],(512-I2[1]))

t = read(FID["/t"])

sct1 = scatter(x=t,y=lEop1)
sct2 = scatter(x=t,y=lEop2)

Plot([sct1,sct2],merge(l_generic,l_compare2))