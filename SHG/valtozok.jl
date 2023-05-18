Nx = 1000
Nt = 600

cry = 4 # GaAs


sigma_t = 1e-12
sigma_x = 5e-3
lambda0 = 10.6e-6
I0 = 60e13

gamma = deg2rad(22)

dz = 1e-6

z_end = 8e-3 + dz


x = range(-sigma_x, sigma_x, Nx) * 6

t = range(-sigma_t, sigma_t, Nt) * 8
