function WannierInterActionBase.Kline(nk::Integer, bse::AbstractBSE; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(nk, bse.TB.lattice, bse.TB.period)
	count(bse.TB.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
function WannierInterActionBase.Kline(nk::Integer, bse::BSESz; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(nk, bse.TB_up.lattice, bse.TB_up.period)
	count(bse.TB.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
function WannierInterActionBase.Kline(point::AbstractVector{<:Pair}, nk::Integer, bse::AbstractBSE; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(point, nk, bse.TB.lattice)
	count(bse.TB.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
function WannierInterActionBase.Kline(point::AbstractVector{<:Pair}, nk::Integer, bse::BSESz; atol::Real = 1e-6, rtol::Real = 1e-5)
	kline = Kline(point, nk, bse.TB_up.lattice)
	count(bse.TB.period) == 3 && _kline_shiftΓ!(kline; atol, rtol)
	return kline
end
