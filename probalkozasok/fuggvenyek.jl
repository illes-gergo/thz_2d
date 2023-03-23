using Symbolics

function deffTHz(cry)
    if cry == 0 # LN
        deff_ = 168e-12
    elseif cry == 2 # ZnTe
        deff_ = 0
    elseif cry == 3 # GaP
        deff_ = 0
    elseif cry == 4 # GaAs
        deff_ = 2 / sqrt(3) * 42.35e-12
    elseif cry == 7 # ZnSe
        deff_ = 0
    end
    return deff_
end

function neo(lambda, T, cry)
    if cry == 4 #GaAs Skauli et al. 2003 0.97-17 um
        l = lambda * 1e6
        a0 = 4.372514
        a = [5.466742 0.0242996 1.957522]
        b = [0.4431307 0.8746453 36.9166]

        n = real(sqrt.(Complex.(a0 .+ 1 .+ a[1] * l .^ 2 ./ (l .^ 2 .- b[1]^2) + a[2] * l .^ 2 ./ (l .^ 2 .- b[2]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 .- b[(3)]^2))))
    end
    return n
end

function ngp(lambda, T, cry)
    lambda1 = lambda * 1e6
    if cry == 4
        @variables l

        a0 = 4.372514
        a = [5.466742 0.0242996 1.957522]
        b = [0.4431307 0.8746453 36.9166]

        n0 = real.(sqrt.(a0 + 1 + a[(1)] * l .^ 2 ./ (l .^ 2 - b[(1)]^2) + a[(2)] * l .^ 2 ./ (l .^ 2 - b[(2)]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 - b[(3)]^2)))

        a = n0 - l * Symbolics.derivative(n0, l)
        #l = lambda1;
        ng = Symbolics.value(substitute(a, l => lambda1))

    end
    return ng
end

function nTHzo(omega, T, cry)
    if cry == 4

        nTHz = real.(sqrt.(er(omega, T, cry)))

    end
    return nTHz
end

function n2value(cry)
    if cry == 4 # GaAs
        n2_ = 5.9e-18
    end
    return n2_
end

function deff(cry)
    if cry == 4 #% GaAs
        deff_ = 65.6e-12
    end
    return deff_
end

function aTHzo(omega, T, cry)
    if cry == 4
        alpha = -2 .* omega / 3e8 .* imag(sqrt.(er(omega, T, cry)))
    end
    return alpha
end
function er(omega, T, cry)
    nu = omega / 2 / pi / 3e8 * 0.01
    if cry == 4 #GaAs
        if T == 300 #ord
            e_inf = 11
            nu_L = 292.1
            nu_T = 268.7
            G = 2.4
            nu = omega / 2 / pi / 3e8 * 1e-2

            er_ = e_inf * (1 .+ (nu_L^2 .- nu_T^2) ./ (nu_T^2 .- nu .^ 2 .+ 1im * G * nu))
        end
    end
    return er_
end