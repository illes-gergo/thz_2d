function createtup(calcdir::String)

    calcs = readdir(calcdir)
    paths = readdir(calcdir, join=true)
    calctuple = ()

    for i in eachindex(calcs)
        calctuple = (calctuple..., (label=calcs[i], value=paths[i]))
    end
    return calctuple
end