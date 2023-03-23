using DSP, Base.Threads, Plots, FourierTools

function conv_manual(a::Array{ComplexF64}, b::Array{ComplexF64})
    conv_res = Array{ComplexF64}(undef, size(a)...)
    @threads for ii = 1:size(a, 1)
        conv_res[ii, :] = sum(a[ii:end, :] .* b[1:end-ii+1, :], dims=1)
    end
    return conv_res
end

function conv_manual(a::Array{Float64}, b::Array{Float64})
    conv_res = Array{Float64}(undef, size(a)...)
    @threads for ii = 1:size(a, 1)
        conv_res[ii, :] = sum(a[ii:end, :] .* b[1:end-ii+1, :], dims=1)
    end
    return conv_res
end



a = vcat(zeros(1000, 1000), rand(1000, 1000) .+ 1im * rand(1000, 1000))
b = vcat(rand(1000, 1000) .+ 1im * rand(1000, 1000), zeros(1000, 1000))

@time A = conv_manual((a), (b))
_, B_ = FourierTools.plan_conv(circshift(a, (1, 0)), reverse((b), dims=(1)), 1)
B_2 = plan_fft(reverse((b), dims=(1)), 1)

@time B = B_(circshift(a, (1, 0)), B_2*(reverse((b), dims=(1))))

plotly()
plot(A[:, 5])
plot!(B[:, 5])