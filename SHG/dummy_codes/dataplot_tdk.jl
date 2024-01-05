using HDF5, PlotlyJS, Interpolations, Base.Threads

e0 = 8.8541878128e-12
c0 = 3e8

include("fuggvenyek.jl")
include("plotfuncs.jl")
include("loutTDK.jl")

calcdir = "/home/illesg/cst/2d-calculations"
calcdir1D = "/home/illesg/cst/calculations/1.0-ps"


config = PlotConfig(
    toImageButtonOptions=attr(
        format="svg", # one of png, svg, jpeg, webp
        filename="custom_image",
        height=480,
        width=640,
        scale=1 # Multiply title/legend/axis/canvas sizes by this factor
    ).fields
)

################ Eltérések ábrázolása ################

#= FID1 = h5open(calcdir * "/20gw1psfull8mm.hdf5", "r")
effic1 = geteffic(FID1)
s1 = scatter(x=(0:400)*20e-3,y=effic1, name="20GW/cm^2")

FID2 = h5open(calcdir * "/40gw1psfull8mm.hdf5", "r")
effic2 = geteffic(FID2)
s2 = scatter(x=(0:400)*20e-3,y=effic2,name="40GW/cm^2")

FID3 = h5open(calcdir * "/60gw1psfull8mm.hdf5", "r")
effic3 = geteffic(FID3)
s3 = scatter(x=(0:400)*20e-3,y=effic3,name="60GW/cm^2")

FID4 = h5open(calcdir * "/80gw1psfull8mm.hdf5", "r")
effic4 = geteffic(FID4)
s4 = scatter(x=(0:400)*20e-3,y=effic4,name="80GW/cm^2")

FID5 = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
effic5 = geteffic(FID5)
s5 = scatter(x=(0:400)*20e-3,y=effic5,name="100GW/cm^2",line_shape="spline")

p = Plot([s1,s2,s3,s4,s5],merge(lout_general,lout_effic))

display(plot(p)) =#

################ 1D ################

#= FID1 = h5open(calcdir1D * "/DB-1ps-20GW")
z = dropdims(read(FID1["/z"])',dims=2)*1e3
effic1 = dropdims(read(FID1["/effic"])',dims=2)*1e2
s1 = scatter(x=z,y=effic1,name="20 GW/cm^2")

FID2 = h5open(calcdir1D * "/DB-1ps-40GW")
effic2 = dropdims(read(FID2["/effic"])',dims=2)*1e2
s2 = scatter(x=z,y=effic2,name="40 GW/cm^2")

FID3 = h5open(calcdir1D * "/DB-1ps-60GW")
effic3 = dropdims(read(FID3["/effic"])',dims=2)*1e2
s3 = scatter(x=z,y=effic3,name="60 GW/cm^2")

FID4 = h5open(calcdir1D * "/DB-1ps-80GW")
effic4 = dropdims(read(FID4["/effic"])',dims=2)*1e2
s4 = scatter(x=z,y=effic4,name="80 GW/cm^2")

FID5 = h5open(calcdir1D * "/DB-1ps-100GW")
effic5 = dropdims(read(FID5["/effic"])',dims=2)*1e2
s5 = scatter(x=z,y=effic5,name="100 GW/cm^2")

p = Plot([s1,s2,s3,s4,s5],merge(lout_general,lout_effic))
 =#
#display(plot(p))


################ Field plots ################

#= FID = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
t = read(FID["/t"]) * 1e12
x = read(FID["/x"]) * 1e3

ETHz = (read(FID["201/ATHz_xt"])) * 2 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) / 1e5

itp = interpolate(abs.(ETHz), BSpline(Cubic()))

iETHz = itp(1:0.5:1024, 1:0.5:2048)

_, I = findmax(abs.(iETHz))
M, ML = extrema(real.(iETHz))

itp = interpolate(real.(ETHz), BSpline(Cubic()))
iETHz = itp(1:0.5:1024, 1:0.5:2048)

p = Plot(contour(x=range(extrema(t)..., 2048), y=range(extrema(x)..., 4096), z=real(circshift(iETHz, (1024 - I[1], 0)))', colorscale=[[0, "#000050"], [0.25, "#0000EF"], [0.5, "#FFFFFF"], [0.75, "#EF0000"], [1, "#500000"]], line_width=0, contours_start=-3300, contours_end=3300, contours_size=33, line_smoothing=1, colorbar=colorbar_field), merge(lout_general, lout_thz))
#display(plot(p))
savefig(p, "100gwfield4mm_new.svg", width=640, height=480) =#

################ Spectrum plots ################

#= FID = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
nu = read(FID["/omega"]) / 1e12 / 2 / pi
nu .-= nu[1]
x = read(FID["/x"]) * 1e3

ATHz = (read(FID["201/ATHz_xo"]))# * 2 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) / 1e5

itp = interpolate(abs.(ATHz), BSpline(Cubic()))

iATHz = itp(1:0.5:1024, 1:0.5:2048)

MX, I = findmax(abs.(iATHz))
M, ML = extrema(abs.(iATHz))

p = Plot(contour(x=range(extrema(nu)..., 2048), y=range(extrema(x)..., 4096), z=abs.(iATHz)' ./ MX, colorscale="Jet", line_width=0, contours_start=0, contours_end=1, contours_size=5e-3, line_smoothing=1, colorbar=colorbar_spectr), merge(lout_general, lout_thz_spectr))
#display(plot(p))
savefig(p, "100gwspect4mm.svg", width=640, height=480) =#

################ Slices ################

#= FID = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
t = read(FID["/t"]) * 1e12
x = read(FID["/x"]) * 1e3

ix = range(extrema(t)..., 20480)
iy = range(extrema(x)..., 4096)

diy = iy[2] - iy[1]

oset_1 = 0.75
oset_2 = -0.75
I_offset_1 = round(Int, oset_1 / diy)
I_offset_2 = round(Int, oset_2 / diy)
println("Residual = $(oset_1-I_offset_1*diy)")

ETHz = (read(FID["201/ATHz_xt"])) * 2 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) / 1e5

itp = interpolate(abs.(ETHz), BSpline(Cubic(Line(OnCell()))))

iETHz = itp(1:0.5:1024, 1:0.5:2048)

_, I = findmax(abs.(iETHz))
M, ML = extrema(real.(iETHz))

itp = interpolate(real(ETHz), BSpline(Cubic(Line(OnCell()))))

iETHz = itp(1:0.5:1024, 1:0.5:2048)

data_1 = real(circshift(iETHz[:, I[2]], (1024 - I[1])))
data_2 = real(circshift(iETHz[:, I[2]+I_offset_1], (1024 - I[1])))
data_3 = real(circshift(iETHz[:, I[2]+I_offset_2], (1024 - I[1])))

idata_1 = interpolate(data_1,BSpline(Cubic(Line(OnCell()))))
idata_2 = interpolate(data_2,BSpline(Cubic(Line(OnCell()))))
idata_3 = interpolate(data_3,BSpline(Cubic(Line(OnCell()))))

s_1 = scatter(x=ix, y=idata_1(collect(1:0.1:2047)), line_shape="spline", name="0 mm")

s_2 = scatter(x=ix, y=idata_2(collect(1:0.1:2047)), line_shape="spline", name="-0.75 mm")

s_3 = scatter(x=ix, y=idata_3(collect(1:0.1:2047)), line_shape="spline", name="0.75 mm")


p = Plot([s_1, s_2, s_3], merge(lout_general, lout_thz1D))

#display(plot(p))
savefig(plot(p), "4mm-slice.svg") =#

################ THz beam intensity profile ################

#= FID = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
x = read(FID["/x"]) * 1e3
z = getSavePoints(FID) * 1e3
thzprofile = getIntensityProfile(FID)
MX, _ = findmax(thzprofile,dims=2)

p = Plot(contour(x=z, y=x, z=thzprofile'./MX', colorscale="Jet", line_width=0, contours_start=0, contours_end=1, contours_size=5e-3, line_smoothing=1, colorbar=colorbar_dist), merge(lout_general,lout_dist))
#display(plot(p))
savefig(p, "100gwintprofilenorm.svg", width=640, height=480) =#

################ THz intensitly profile slice ################

#= FID = h5open(calcdir * "/100gw1psfull8mm.hdf5", "r")
x = read(FID["/x"]) * 1e3
z = getSavePoints(FID) * 1e3
thzprofile = getIntensityProfile(FID)
#MX, _ = findmax(thzprofile, dims=2)
thzSlice = thzprofile[291, :]
MX, _ = findmax(thzSlice)

p = Plot(scatter(x=x, y=thzSlice/MX, line_shape="spline"), merge(lout_general, lout_int1D))
#display(plot(p))
savefig(p, "100gw6mmintprofile.svg", width=640, height=480) =#

################ THz Centroid Frequency Profile ################

FID = h5open(calcdir * "/20gw1psfull8mm.hdf5", "r")
x = read(FID["/x"]) * 1e3
z = getSavePoints(FID) * 1e3
thzprofile = getCentralFrequency(FID)
#MX, _ = findmax(thzprofile, dims=2)
#thzSlice = thzprofile[291, :]
MX, _ = findmax(thzprofile)

p = Plot(contour(x=z, y=x, z=thzprofile', colorscale="Jet", line_width=0, line_smoothing=1, contours_start=0,contours_end=3,contours_size=0.02, colorbar=colorbar_dist_spect), merge(lout_general, lout_dist_spect))

#display(plot(p))
savefig(p, "dist_spect_20gw.png", width=640, height=480, scale = 5)