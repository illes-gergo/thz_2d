Nx = 2048 #2048
Nt = 1024 # 1024

cry = 4 # GaAs


sigma_t = 1e-12
sigma_x = 2e-3
lambda0 = 10.6e-6
I0 = 1e13

STR = "próbaszámolás-koherencia-ss"

#gamma = deg2rad(22)
gamma = acos(ngp(lambda0, 300, cry) / nTHzo(1.5e12 * 2 * pi, 300, cry))

dz = 1e-6

z_end = 0.5e-3 + dz


x = range(-sigma_x, sigma_x, Nx) * 5

t = range(-sigma_t, sigma_t, Nt) * 20