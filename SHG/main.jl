using LazyGrids, FFTW, FourierTools, Base.Threads, Plots, Dates, HDF5

# FFT -> /omegaMAX ; IFFT -> * omegaMAX

default(levels=100)
gr()
include("gauss_impulzus.jl")
include("diffegy_megoldo.jl")
include("differencial_egyenletek.jl")
include("fuggvenyek.jl")
include("valtozok.jl")

const c0 = 3e8
khi_eff = 2 * deffTHz(cry)
d_eff = deff(cry)
const e0 = 8.854187817e-12

SHG_SHIFT = floor(Int, Nt / 4)

tMax = t[end] - t[1]
dt = t[2] - t[1]

xMax = x[end] - x[1]
dx = x[2] - x[1]

dOmega = 2 * pi / tMax
omega = range(0, length=Nt, step=dOmega)

omegaMax = omega[end] - omega[1]
omega = omega .- omegaMax ./ 2

dkx = 2 * pi / xMax
kx = range(0, length=Nx, step=dkx) .- Nx / 2 * dkx
kx0 = sin(gamma) ./ lambda0 .* 2 .* pi * neo(lambda0, 300, cry)

kxMax = kx[end] - kx[1]

omega_diff = round(Int, 2 * pi * c0 / lambda0 / dOmega)
omega0 = omega_diff .* dOmega
lambda0 = 2 * pi * c0 / omega0


(ct, cx) = ndgrid(t, x)
(comega_, ckx) = ndgrid(omega, kx)

comega = comega_ .+ omega0

comegaSHG = comega_ .+ 2 * omega0

comegaTHz = comega_ .- omega[1]

clambda = c0 ./ comega * 2 * pi

cLambdaSHG = c0 ./ comegaSHG * 2 * pi

n = neo(clambda, 300, cry)

k_omega = n .* comega ./ c0
kx_omega = real.(k_omega .* sin(gamma))
kz_omega = real.(k_omega .* cos(gamma))

k_omegaSHG = neo(cLambdaSHG, 300, cry) .* comegaSHG ./ c0
kx_omegaSHG = k_omegaSHG .* sin(gamma)
kz_omegaSHG = k_omegaSHG .* cos(gamma)

nTHz = nTHzo(comegaTHz, 300, cry)

k_omegaTHz = nTHz .* comegaTHz ./ c0
kz_omegaTHz = real.(sqrt.(Complex.(k_omegaTHz .^ 2 - ckx .^ 2)))

E0 = sqrt(2 * I0 / neo(lambda0, 300, cry) / e0 / c0)

Axt = gauss_impulzus_omega0(E0, sigma_t, sigma_x, lambda0, gamma, ct, cx)

fft_t_o = plan_fft(Axt, 1)
fft_x_kx = plan_fft(Axt, 2)
ifft_o_t = plan_ifft(Axt, 1)
ifft_kx_x = plan_ifft(Axt, 2)

padding = zeros(Nt, Nx)
plan_fast_conv(Axt, Axt)

alpha = aTHzo(comegaTHz, 300, cry)

Axo = fftshift(fft_t_o * Axt, 1) ./ omegaMax .* exp.(+1im .* kx_omega .* cx)
Akxo = fftshift(fft_x_kx * Axo / kxMax, 2)
z = Array{Float64}(undef, floor(Int, z_end / dz))

ATHz_kx_o = zeros(size(Akxo))

ASH = copy(ATHz_kx_o)

A_kompozit = cat(Akxo, ATHz_kx_o, ASH, dims=3)

z[1] = 0

global plotInteraction::Bool = false
#STR = Dates.format(now(), "YYYY-MM-DD hh-mm-ss")
STR = "próbaszámolás.hdf5"
global Axo_prew = zeros(size(Axo))
FID = h5open(STR * ".hdf5", "w")
entryCounter::Int = 1;
#STR = "elojel_minusz"
for ii in 1:(length(z)-1)
    global A_kompozit, z[ii+1] = RK4M(thz_feedback_n2_SHG, z[ii], A_kompozit, dz)
    if mod(ii, 11) == 0
        global plotInteraction = true
    else
        global plotInteraction = false
    end

    #if (mod(ii, 100) == 0 || ii == 1 ) && false
    if ii == length(z) - 1 || mod(ii, 20) == 0 || ii == 1

        global Aop_kx_o = A_kompozit[:, :, 1]
        global Axo = ifft_kx_x * ifftshift(Aop_kx_o, 2) .* kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* z[ii+1])
        global Axt = ifft_o_t * ifftshift(Axo .* omegaMax, 1)
        global Aop_kx_oSH = A_kompozit[:, :, 3]
        global AxoSH = ifft_kx_x * ifftshift(Aop_kx_oSH, 2) .* kxMax .* exp.(-1im .* kx_omegaSHG .* cx - 1im .* kz_omegaSHG .* z[ii+1])
        global AxtSH = ifft_o_t * ifftshift(AxoSH .* omegaMax, 1)

        p1 = heatmap(x, t, abs.(Axt .* exp.(1im .* omega0 .* ct) .* 1e-8), linewidth=0, colormap=:jet)
        p3 = heatmap(x, omega, abs.(Axo - Axo_prew), colormap=:jet)
        p4 = heatmap(x, t, abs.(AxtSH .* exp.(2im .* omega0 .* ct) .* 1e-8), linewidth=0, colormap=:jet)
        global ATHz_kx_o = A_kompozit[:, :, 2]
        global ATHz_xo = ifft_kx_x * ifftshift(ATHz_kx_o .* exp.(-1im .* k_omegaTHz .* z[ii+1]), 2) .* kxMax
        global ATHz_xt = ifft_o_t * ATHz_xo * omegaMax
        p2 = heatmap(x, t, real.(ATHz_xt) * 1e-5, linewidth=0, colormap=:jet)
        global _, max_indices = findmax(abs.(Axt))
        (scatter!([x[max_indices[2]]], [t[max_indices[1]]]))
        display(plot(p1, p2, p3, p4, layout=(2, 2), size=[1200, 900]))
        global Axo_prew = copy(Axo)
        #display(heatmap(x, t, abs.(ATHz_kx_o), linewidth=0, colormap=:jet))
        FID["/"*string(entryCounter)*"/Eop"] = Axt
        FID["/"*string(entryCounter)*"/Aop"] = Axo
        FID["/"*string(entryCounter)*"/ASH"] = AxoSH
        FID["/"*string(entryCounter)*"/ATHz_xo"] = ATHz_xo
        FID["/"*string(entryCounter)*"/ATHz_xt"] = ATHz_xt
        global entryCounter += 1
    end
    display(ii)
end
FID["/maxEntry"] = entryCounter - 1
FID["/gamma"] = rad2deg(gamma)
FID["/z"] = z
FID["/omega"] = omega
FID["/omega0"] = omega0
FID["/x"] = collect(x)
FID["/t"] = collect(t)
close(FID)