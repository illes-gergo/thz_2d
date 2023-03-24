function imp_terjedes(t, Y)
    dAdz = -1im .* kx_omega .* ckx ./ kz_omega .* Y + 1im .* ckx .^ 2 ./ 2 ./ kz_omega .* Y
    return dAdz
end

function thz_generation(t, Y)
    Eop = ifft_kx_x * ifftshift(Y, 2) * kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* t)
    conv_part = fast_forward_convolution(Eop, conj(Eop)) * e0 * d_eff * dOmega
    return fftshift(fft_x_kx * (conv_part), 2) ./ kxMax
end

function thz_egyszeru(t, Y)
    Aop = Y[:, :, 1]
    ATHz = Y[:, :, 2]

    dAopdz = @spawn imp_terjedes(t, Aop)
    dTHz_gen = @spawn begin
        temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* k_omegaTHz .* t) - alpha / 2 .* ATHz
        temp_val[isnan.(temp_val)] .= 0
        return temp_val
    end

    return cat(fetch(dAopdz), fetch(dTHz_gen), dims=3)
end

function plan_fast_conv(a, b)
    a_ = vcat(padding, a)
    b_ = vcat(b, padding)
    global _, fast_conv_plan = FourierTools.plan_conv(a_, b_, 1)
    global fast_conv_fft_plan = plan_fft(b_, 1)
end

function fast_forward_convolution(a, b)
    a_ = vcat(padding, a)
    b_ = vcat(b, padding)
    return fast_conv_plan(circshift(a_, (1, 0)), fast_conv_fft_plan * (reverse((b_), dims=(1))))[floor(Int, end / 2)+1:end, :]
end

function fast_backward_convolution(a, b)
    a_ = vcat(padding, a)
    b_ = vcat(b, padding)
    return reverse(fast_conv_plan(circshift(reverse(a_, dims=1), (1, 0)), fast_conv_fft_plan * (reverse((b_), dims=(1)))), dims=1)[floor(Int, end / 2)+1:end, :]
end