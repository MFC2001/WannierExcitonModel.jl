
function ExtractU(qpoints::AbstractVector{<:ReducedCoordinates}, band₁::AbstractVector{<:Eigen}, band₂::AbstractVector{<:Eigen}, bse::BSEspinless; η)
	length(qpoints) == length(band₁) == length(band₂) || error("Mismatched qpoints and band.")
	ijRvck = _uijR_ψvck(bse, η)
	band₁_u = similar(band₁)
	band₂_u = similar(band₂)
	for (qi, q) in enumerate(qpoints)
		q = _BSE_preprocess_eleband_q!(bse, q, Val(bse.isqgrid))
		BM = ijRvck(bse.bandk, bse.bandkq, q)
		band₁_u[qi] = Eigen(copy(band₁[qi].values), BM * band₁[qi].vectors)
		band₂_u[qi] = Eigen(copy(band₂[qi].values), BM * band₂[qi].vectors)
	end
	return band₁_u, band₂_u
end

function ExtractU(q::ReducedCoordinates, band₁::Eigen, band₂::Eigen, bse::BSEspinless; η)
	ijRvck = _uijR_ψvck(bse, η)
	q = _BSE_preprocess_eleband_q!(bse, q, Val(bse.isqgrid))
	BM = ijRvck(bse.bandk, bse.bandkq, q)
	band₁_u = Eigen(copy(band₁.values), BM * band₁.vectors)
	band₂_u = Eigen(copy(band₂.values), BM * band₂.vectors)
	return band₁_u, band₂_u
end

