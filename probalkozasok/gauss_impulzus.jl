function gauss_impulzus(E0, sigma_t, sigma_x, omega0, lambda0, gamma, t, x)
    c0 = 3e8
    return E0 .* exp.(-2 .* log(2) .* t .^ 2 ./ sigma_t .^ 2) .*
           exp.(-x .^ 2 ./ sigma_x .^ 2) .* exp.(1im .* omega0 .* t) .* exp.(-1im .* sin(gamma) .* x ./ lambda0 .* 2 .* pi * neo(lambda0, 300, cry))
end
function gauss_impulzus_omega0(E0, sigma_t, sigma_x, lambda0, gamma, t, x)
    c0 = 3e8
    return E0 .* exp.(-2 .* log(2) .* t .^ 2 ./ sigma_t .^ 2) .*
           exp.(-x .^ 2 ./ sigma_x .^ 2) .* exp.(-1im .* sin(gamma) .* x ./ lambda0 .* 2 .* pi * neo(lambda0, 300, cry))
end