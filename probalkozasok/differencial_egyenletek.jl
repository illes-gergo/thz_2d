function imp_terjedes(t, Y)
    dAdz = -1im .* kx_omega .* ckx ./ kz_omega .* Y + 1im .* ckx.^2 ./ 2 ./ kz_omega .* Y
    return dAdz
end