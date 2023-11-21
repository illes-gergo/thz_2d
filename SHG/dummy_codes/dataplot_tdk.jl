using HDF5, PlotlyJS

e0 = 8.8541878128e-12
c0 = 3e8

include("fuggvenyek.jl")
include("plotfuncs.jl")
include("loutTDK.jl")

calcdir = "/home/illesg/cst/2d-calculations"
calcdir1D = "/home/illesg/cst/calculations/2.0-ps"

#= FID1 = h5open(calcdir * "/20gw2psfull8mm.hdf5", "r")
effic1 = geteffic(FID1)
s1 = scatter(x=(0:400)*20e-3,y=effic1, name="20GW/cm^2")

FID2 = h5open(calcdir * "/40gw2psfull8mm.hdf5", "r")
effic2 = geteffic(FID2)
s2 = scatter(x=(0:400)*20e-3,y=effic2,name="40GW/cm^2")

FID3 = h5open(calcdir * "/60gw2psfull8mm.hdf5", "r")
effic3 = geteffic(FID3)
s3 = scatter(x=(0:400)*20e-3,y=effic3,name="60GW/cm^2")

FID4 = h5open(calcdir * "/80gw2psfull8mm.hdf5", "r")
effic4 = geteffic(FID4)
s4 = scatter(x=(0:400)*20e-3,y=effic4,name="80GW/cm^2")

FID5 = h5open(calcdir * "/100gw2psfull8mm.hdf5", "r")
effic5 = geteffic(FID5)
s5 = scatter(x=(0:400)*20e-3,y=effic5,name="100GW/cm^2")

p = Plot([s1,s2,s3,s4,s5],merge(lout_general,lout_effic))

display(plot(p)) =#

# 1D

FID1 = h5open(calcdir1D * "/DB-2ps-20GW")
z = dropdims(read(FID1["/z"])',dims=2)*1e3
effic1 = dropdims(read(FID1["/effic"])',dims=2)*1e2
s1 = scatter(x=z,y=effic1,name="20 GW/cm^2")

FID2 = h5open(calcdir1D * "/DB-2ps-40GW")
effic2 = dropdims(read(FID2["/effic"])',dims=2)*1e2
s2 = scatter(x=z,y=effic2,name="40 GW/cm^2")

FID3 = h5open(calcdir1D * "/DB-2ps-60GW")
effic3 = dropdims(read(FID3["/effic"])',dims=2)*1e2
s3 = scatter(x=z,y=effic3,name="60 GW/cm^2")

FID4 = h5open(calcdir1D * "/DB-2ps-80GW")
effic4 = dropdims(read(FID4["/effic"])',dims=2)*1e2
s4 = scatter(x=z,y=effic4,name="80 GW/cm^2")

FID5 = h5open(calcdir1D * "/DB-2ps-100GW")
effic5 = dropdims(read(FID5["/effic"])',dims=2)*1e2
s5 = scatter(x=z,y=effic5,name="100 GW/cm^2")

p = Plot([s1,s2,s3,s4,s5],merge(lout_general,lout_effic))

display(plot(p))