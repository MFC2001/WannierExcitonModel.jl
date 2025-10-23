
include("./SU2.jl")
# include("./SP.jl")
# include("./wannier_pp.jl")

include("./mmn.jl")
include("./amn.jl")
include("./guess.jl")

export BSEwannier

function BSEwannier(qgrid::MonkhorstPack, bse::AbstractBSE; kwargs...)
	return BSEwannier(RedKgrid(qgrid), bse; kwargs...)
end
