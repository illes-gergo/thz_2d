Nx = 500
Nt = 500

cry = 0 # GaAs


sigma_t = 0.5e-12
sigma_x = 2e-3
lambda0 = 1030e-9
I0 = 100e13

gamma = deg2rad(63.5)

dz = 5e-7

z_end = 1e-3


x = range(-sigma_x, sigma_x, Nx) * 8

t = range(-sigma_t, sigma_t, Nt) * 8
