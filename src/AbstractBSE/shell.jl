function WannierInterActionBase.Kline(nk::Integer, bse::AbstractBSE; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(nk, bse.lattice, bse.period)
	count(bse.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
function WannierInterActionBase.Kline(point::AbstractVector{<:Pair}, nk::Integer, bse::AbstractBSE; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(point, nk, bse.lattice)
	count(bse.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
