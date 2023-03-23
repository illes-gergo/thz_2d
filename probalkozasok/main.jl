using LazyGrids, FFTW, Plots, FourierTools, Base.Threads

# FFT -> /omegaMAX ; IFFT -> * omegaMAX

default(size=(800, 800), levels=50)

include("valtozok.jl")
include("gauss_impulzus.jl")
include("diffegy_megoldo.jl")
include("differencial_egyenletek.jl")
include("fuggvenyek.jl")

const c0 = 3e8
d_eff = deffTHz(cry)
e0 = 8.854187817e-12

tMax = t[end] - t[1]
dt = t[2] - t[1]

xMax = x[end] - x[1]
dx = x[2] - x[1]

dOmega = 2 * pi / tMax
omega = range(0, length=N, step=dOmega)

omegaMax = omega[end] - omega[1]
omega = omega .- omegaMax ./ 2

dkx = 2 * pi / xMax
kx = range(0, length=N, step=dkx) .- N / 2 * dkx

kxMax = kx[end] - kx[1]

omega_diff = round(Int, 2 * pi * c0 / lambda0 / dOmega)
omega0 = omega_diff .* dOmega
lambda0 = 2 * pi * c0 / omega0


(ct, cx) = ndgrid(t, x)
(comega_, ckx) = ndgrid(omega, kx)

comega = comega_ .+ omega0

comegaTHz = comega_ .- omega[1]

clambda = c0 ./ comega * 2 * pi


#= neo(lambda::Array, T, cry) = ones(size(lambda)...)
neo(lambda::Number, T, cry) = 1
 =#
n = neo(clambda, 300, cry)

k_omega = n .* comega ./ c0
k_omega[1, :] = k_omega[2, :]
kx_omega = k_omega .* sin(gamma)
kz_omega = k_omega .* cos(gamma)

nTHz = nTHzo(comegaTHz, 300, cry)

k_omegaTHz = n .* comegaTHz ./ c0
kz_omegaTHz = sqrt.(Complex.(k_omegaTHz .^ 2 - ckx .^ 2))
kz_omegaTHz -= 1im * imag(kz_omegaTHz)

Axt = gauss_impulzus_omega0(E0, sigma_t, sigma_x, lambda0, gamma, ct, cx)

#= display(contourf(x, t, abs.(Axt), colormap=:jet, linewidth=0))
display(contourf(x, t, real.(Axt), colormap=:jet, linewidth=0)) =#


fft_t_o = plan_fft(Axt, 1)
fft_x_kx = plan_fft(Axt, 2)
ifft_o_t = plan_ifft(Axt, 1)
ifft_kx_x = plan_ifft(Axt, 2)

padding = zeros(size(Axt))
plan_fast_conv(Axt, Axt)

Axo = fftshift(fft_t_o * Axt, 1) ./ omegaMax .* exp.(+1im .* kx_omega .* cx)
Akxo = fftshift(fft_x_kx * Axo / kxMax, 2)
#= Axo2 = ifft_kx_x * ifftshift(Akxo .* kxMax, 2) .* exp.(-1im .* kx_omega .* cx)
Axt2 = ifft_o_t * ifftshift(Axo2, 1) .* omegaMax
display(contourf(x, omega, abs.(Axo), colormap=:jet, linewidth=0))
display(plot(omega, abs.(Axo[:, Int(end / 2)])))
display(contourf(kx, omega, abs.(Akxo), colormap=:jet, linewidth=0, xlim=[-kxMax, kxMax] / 10))
display(plot(omega, abs.(Akxo[:, Int(end / 2)])))
display(contourf(x, t, real.(Axt2 .* exp.(1im .* omega0 .* t)), colormap=:jet, linewidth=0))
 =#
z = Array{Float64}(undef, N * 50)

ATHz_kx_o = zeros(size(Akxo))

A_kompozit = cat(Akxo, ATHz_kx_o, dims=3)

z[1] = 0;


#= let Akxo = Akxo
    for ii in 1:(length(z)-1)
        Akxo, z[ii+1] = RK4M(imp_terjedes, z[ii], Akxo, dz)
        if mod(ii, 500) == 0
            #display(contourf(kx, omega, abs.(Akxo), linewidth=0, xlim=[-kxMax, kxMax] / 2, colormap=:jet))
            local Axo = ifft_kx_x * ifftshift(Akxo, 2) .* kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* z[ii+1])
            local Axt = ifft_o_t * ifftshift(Axo .* omegaMax, 1)
            display(contourf(x, t, real.(Axt .* exp.(1im .* omega0 .* t)), linewidth=0, colormap=:jet))
        end
        display(ii)
    end
end =#

let A_kompozit = A_kompozit
    for ii in 1:(length(z)-1)
        A_kompozit, z[ii+1] = RK4M(thz_egyszeru, z[ii], A_kompozit, dz)
        if mod(ii, 100) == 0 || ii == 1
            local Aop_kx_o = A_kompozit[:, :, 1]
            #display(contourf(kx, omega, abs.(Akxo), linewidth=0, xlim=[-kxMax, kxMax] / 2, colormap=:jet))
            local Axo = ifft_kx_x * ifftshift(Aop_kx_o, 2) .* kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* z[ii+1])
            local Axt = ifft_o_t * ifftshift(Axo .* omegaMax, 1)
            display(contourf(x, t, real.(Axt .* exp.(1im .* omega0 .* t)), linewidth=0, colormap=:jet))
            local ATHz_kx_o = A_kompozit[:, :, 2]
            local ATHz_xo = ifft_kx_x * ifftshift(ATHz_kx_o, 2) .* kxMax .* exp.(-1im .* kz_omegaTHz .* z[ii+1])
            local ATHz_xt = ifft_o_t * ATHz_xo * omegaMax
            contourf(x, t, real.(ATHz_xt), linewidth=0, colormap=:jet)
            local _, max_indices = findmax(abs.(Axt))
            display(scatter!([x[max_indices[2]]], [t[max_indices[1]]]))
            #display(contourf(x, t, abs.(ATHz_kx_o), linewidth=0, colormap=:jet))
        end
        display(ii)
    end
end