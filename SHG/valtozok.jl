Nx = 2000
Nt = 800

cry = 4 # GaAs


sigma_t = 1e-12
sigma_x = 5e-3
lambda0 = 10.6e-6
I0 = 60e13

#gamma = deg2rad(22)
gamma = acos(ngp(lambda0, 300, cry) / nTHzo(0.5e12 * 2 * pi, 300, cry))

dz = 1e-6

z_end = 8e-3 + dz


x = range(-sigma_x, sigma_x, Nx) * 10

t = range(-sigma_t, sigma_t, Nt) * 20
