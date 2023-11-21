function I_opMax(inds)
    out = zeros(size(inds))
    for i in eachindex(inds)
        out[i], _ = findmax(abs.(read(FID["/$(i)/Eop"])) .^ 2 * e0 * c0 / 2 * neo(10.6e-6,300,4) / 1e13)
    end
    return out
end

function geteffic(FID)
    En0 = sum(read(FID["/En"])[1, :])
    EnTHz = sum(read(FID["/EnTHz"]), dims=2)

effic = EnTHz ./ En0

effic = dropdims(effic, dims=2) * 100 * 4 * nTHzo(1.5e12 * 2 * pi, 300, 4) / (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) .^ 2
return effic

end