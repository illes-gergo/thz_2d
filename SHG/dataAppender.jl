using HDF5, Base.Threads, DSP
include("fuggvenyek.jl")

calcdir = "/home/illesg/cst/2d-calculations/"


#readdir(calcdir)
dirs = ["20gw2psfull8mm-narrow-alt",
	"40gw2psfull8mm-narrow-alt",
	"60gw2psfull8mm-narrow-alt",
	"80gw2psfull8mm-narrow-alt",
	"100gw2psfull8mm-narrow-alt"].*".hdf5"
for dir in dirs
    
    FID = h5open(calcdir * dir,"r+")

    if haskey(FID,"En") delete_object(FID,"En") end
    if haskey(FID,"EnTHz") delete_object(FID,"EnTHz") end
    if haskey(FID,"EnSH") delete_object(FID,"EnSH") end

    En = zeros(401,2048)
    EnTHz = zeros(401,2048)
    EnSH = zeros(401,2048)
    En0 = sum(abs.(read(FID["1/Aop"])).^2,dims=1) * neo(10.6e-6,300,4)
    for i in 1:401
        EnTHz[i,:] = sum(abs.(read(FID[string(i)*"/ATHz_xo"])).^2,dims=1) .* nTHzo(1.5e12*2*pi,300,4)
        EnSH[i,:] = sum(abs.(read(FID[string(i)*"/ASH"])).^2,dims=1) .* neo(10.6e-6/2,300,4)
        En[i,:] = sum(abs.(read(FID[string(i)*"/Aop"])).^2,dims=1) * neo(10.6e-6,300,4)
        println("ITER $(i)")
    end
    FID["En"] = En
    FID["EnTHz"] = EnTHz
    FID["EnSH"] = EnSH
    println("DB $(dir)")
end
