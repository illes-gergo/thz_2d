function I_opMax(inds)
    out = zeros(size(inds))
    for i in eachindex(inds)
        out[i], _ = findmax(abs.(read(FID["/$(i)/Eop"])) .^ 2 * e0 * c0 / 2 * neo(10.6e-6, 300, 4) / 1e13)
    end
    return out
end

function geteffic(FID)
    En0 = sum(read(FID["/En"])[1, :])
    EnTHz = sum(read(FID["/EnTHz"]), dims=2)

    effic = EnTHz ./ En0

    effic = dropdims(effic, dims=2) * 100 * 4 * nTHzo(1.5e12 * 2 * pi, 300, 4) ./ (1 + nTHzo(1.5e12 * 2 * pi, 300, 4)) .^ 2
    return effic

end

function getSavePoints(FID)
    z_in_DB = read(FID["z"])
    return range(0, z_in_DB[end], read(FID["/maxEntry"]))
end

function getIntensityProfile(FID)
    return read(FID["EnTHz"])
end

function getCentralFrequency(FID)
    omega = read(FID["omega"])
    nu = (omega .- omega[1]) / 2 / pi / 1e12
    maxentry = read(FID["maxEntry"])
    retval = zeros(maxentry, 2048)
    for i in 1:maxentry
        spectr = abs.(read(FID["$i/ATHz_xo"])) .^ 2
        @threads for j in 1:2048
            retval[i, j] = sum(nu .* spectr[:, j]) / sum(spectr[:, j])
        end
    end
    return retval
end