
include("./spinless.jl")
include("./spinful.jl")

function ExtractU(qpoints::AbstractVector{<:ReducedCoordinates}, band::AbstractVector{<:Eigen}, bse::AbstractBSE; η)
	length(qpoints) == length(band) || error("Mismatched qpoints and band.")
	ijRvck = _uijR_ψvck(bse, η)
	band_u = similar(band)
	for (qi, q) in enumerate(qpoints)
		q = _BSE_preprocess_eleband_q!(bse, q, Val(bse.isqgrid))
		BM = ijRvck(bse.bandk, bse.bandkq, q)
		band_u[qi] = Eigen(copy(band[qi].values), BM * band[qi].vectors)
	end
	return band_u
end

function ExtractU(q::ReducedCoordinates, band::Eigen, bse::AbstractBSE; η)
	ijRvck = _uijR_ψvck(bse, η)
	q = _BSE_preprocess_eleband_q!(bse, q, Val(bse.isqgrid))
	BM = ijRvck(bse.bandk, bse.bandkq, q)
	band_u = Eigen(copy(band.values), BM * band.vectors)
	return band_u
end
