
include("./spinless.jl")
# include("./SP.jl")
# include("./wannier_pp.jl")

# include("./mmn.jl")
include("./amn.jl")
include("./guess.jl")

export BSE_wannier

function BSE_wannier(qgrid::MonkhorstPack, bse::AbstractBSE; kwargs...)
	return BSE_wannier(RedKgrid(qgrid), bse; kwargs...)
end
