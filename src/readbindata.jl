function readbindata(filename::String="dat/PT-ZB18a_1-1.bin")
    open(filename, "r") do io
        kount = Int(read(io,Int32))
        buf = Vector{Float64}(undef,kount)
        time = Vector{Float64}(undef,kount); time = read!(io, time)
        eccentricity = Vector{Float64}(undef,kount); eccentricity = read!(io, eccentricity)
        obliquity = Vector{Float64}(undef,kount); obliquity = read!(io, obliquity)
        precession = Vector{Float64}(undef,kount); precession = read!(io, precession)
        lpx = Vector{Float64}(undef,kount); lpx = read!(io, lpx)
        climatic_precession = Vector{Float64}(undef,kount); climatic_precession = read!(io, climatic_precession)

        return DataFrame(
            time = time / 365250., # days to kyr
            eccentricity = eccentricity,
            obliquity = obliquity,
            precession = precession,
            lpx = lpx,
            climatic_precession = climatic_precession
        )
    end
end
