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

function conv_manual_backwards(a::Array{ComplexF64}, b::Array{ComplexF64})
    conv_res = Array{ComplexF64}(undef, size(a)...)
    for ii = 1:size(a, 1)
        conv_res[ii, :] = sum(a[ii:-1:1, :] .* b[1:ii, :], dims=1)
    end
    return conv_res
end

a = rand(10, 10) .+ 1im * rand(10, 10)
b = rand(10, 10) .+ 1im * rand(10, 10)
a_ = vcat(zeros(10, 10), a)
b_ = vcat(b, zeros(10, 10))

# @time A = conv_manual((a), (b))
_, B_ = FourierTools.plan_conv(a_, b_, 1)
B_2 = plan_fft(b_, 1)

#@time B = B_(circshift(a_, (1, 0)), B_2 * (reverse((b_), dims=(1))))[floor(Int, end / 2)+1:end, :]


A = conv_manual_backwards(a, b)
B = reverse(conv_manual(reverse(a, dims=1), (b)), dims=1)
C = reverse(B_(circshift(reverse(a_, dims=1), (1, 0)), B_2 * (reverse((b_), dims=(1)))),dims=1)[floor(Int, end / 2)+1:end, :]
#plotly()
plot(A[:, 5])
#plot!(B[:, 5])
plot!(C[:, 5])