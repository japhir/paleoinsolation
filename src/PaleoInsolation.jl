module PaleoInsolation

using DataFrames

include("insolation.jl")
include("readbindata.jl")

export insolation
export readbindata

end # module
