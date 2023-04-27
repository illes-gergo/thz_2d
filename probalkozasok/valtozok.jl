Nx = 500
Nt = 500

cry = 4 # GaAs


sigma_t = 1e-12
sigma_x = 2e-3
lambda0 = 10.6e-6
I0 = 100e13

gamma = deg2rad(21)

dz = 1e-6

z_end = 1e-3 + dz


x = range(-sigma_x, sigma_x, Nx) * 8

t = range(-sigma_t, sigma_t, Nt) * 8
