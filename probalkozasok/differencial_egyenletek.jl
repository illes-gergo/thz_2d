function imp_terjedes(t, Y)
    dAdz = -1im .* kx_omega .* ckx ./ kz_omega .* Y + 1im .* ckx .^ 2 ./ 2 ./ kz_omega .* Y
    return dAdz
end

function thz_generation(t, Y)
    Eop = ifft_kx_x * ifftshift(Y,2) * kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* t)
end

function thz_egyszeru(t, Y)
    Aop = Y[:, :, 1]
    ATHz = Y[:, :, 2]

    dAopdz = @spawn imp_terjedes(t, Aop)

    return cat(dAopdz.result, dATHz.result, dims=3)
end