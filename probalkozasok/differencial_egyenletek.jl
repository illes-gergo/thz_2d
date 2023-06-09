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
        temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
        temp_val[isnan.(temp_val)] .= 0
        return temp_val
    end
    wait.([dAopdz, dTHz_gen])
    return cat(dAopdz.result, dTHz_gen.result, dims=3)
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

function thz_cascade(t, Aop, ATHz)
    Eop = @spawn ifft_kx_x * ifftshift(Aop, 2) * kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* t)
    ETHz = @spawn ifft_kx_x * ifftshift(ATHz * kxMax .* exp.(-1im .* kz_omegaTHz .* t), 2)
    wait.([Eop, ETHz])
    temp_val1 = @spawn e0 .* d_eff .* fast_forward_convolution(Eop.result, conj(ETHz.result))
    temp_val2 = @spawn e0 .* d_eff .* fast_backward_convolution(Eop.result, ETHz.result)
    wait.([temp_val1, temp_val2])
    return fftshift(fft_x_kx * ((temp_val1.result .+ temp_val2.result) .* exp.(+1im .* kx_omega .* cx)) / kxMax .* dOmega,2)
end

function thz_feedback(t, Y)
    Aop = Y[:, :, 1]
    ATHz = Y[:, :, 2]
    dAop_lin = @spawn imp_terjedes(t, Aop)
    dTHz_gen = @spawn begin
        temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
        temp_val[isnan.(temp_val)] .= 0
        return temp_val
    end
    dAopCsc = @spawn thz_cascade(t, Aop, ATHz) .* exp.(1im .* kz_omega .* t)

    wait.([dAop_lin, dTHz_gen, dAopCsc])
    #= println("Wait ok")
    display(plot(heatmap(abs.(dAopCsc.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2))))
    display(plot(heatmap(abs.(Aop))))
    println("Plot ok")
    #p2 = plot(contour(z=abs.(Aop)))
    readline()
    println("readline() ok") =#
    return cat(dAop_lin.result .- dAopCsc.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2, dTHz_gen.result, dims=3)
end

function n2calc(t, Aop)
    mult1 = @spawn -2 .* kz_omega .* e0 .* c0 .^ 2 ./ comega .^ 2
    Eop = @spawn ifft_o_t * ifftshift(ifft_kx_x * ifftshift(Aop, 2) * kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* t), 1) * omegaMax
    mult2 = e0 * omega0 * neo(lambda0, 300, cry) * n2value(cry) / 2
    wait(Eop)
    aEop2 = (Eop.result) .* conj(Eop.result)
    wait(mult1)
    return fftshift(fft_x_kx * (mult1.result .* fftshift(fft_t_o * (mult2 .* aEop2 .* Eop.result), 1) ./ omegaMax .* exp.(+1im .* kx_omega .* cx) ./ kxMax), 2)
end

function thz_feedback_n2(t, Y)
    Aop = Y[:, :, 1]
    ATHz = Y[:, :, 2]
    dAop_lin = @spawn imp_terjedes(t, Aop)
    dTHz_gen = @spawn begin
        temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
        temp_val[isnan.(temp_val)] .= 0
        return temp_val
    end
    dAopCsc = @spawn thz_cascade(t, Aop, ATHz) #.* exp.(1im .* kz_omega .* t)
    dAopn2 = @spawn n2calc(t, Aop)

    sum_dAop = @spawn begin
        wait.([dAopCsc, dAopn2])
        return (dAopCsc.result .- dAopn2.result) .* exp.(1im .* kz_omega .* t)
    end
    #= println("Wait ok")
    display(plot(heatmap(abs.(dAopCsc.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2))))
    display(plot(heatmap(abs.(Aop))))
    println("Plot ok")
    #p2 = plot(contour(z=abs.(Aop)))
    readline()
    println("readline() ok") =#
    wait.([sum_dAop, dTHz_gen, dAop_lin])
    return cat(dAop_lin.result .- sum_dAop.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2, dTHz_gen.result, dims=3)
end