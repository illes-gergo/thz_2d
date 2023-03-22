function RK4M(f::Function, t, Y, step)
    k1 = f(t, Y)
    k2 = f(t .+ step ./ 2, Y .+ step ./ 2 .* k1)
    k3 = f(t .+ step ./ 2, Y .+ step ./ 2 .* k2)
    k4 = f(t .+ step, Y .+ step .* k3)
    Y += step .* (k1 .+ 2 * k2 .+ 2 .* k3 .+ k4) ./6
    return Y, t + step
end