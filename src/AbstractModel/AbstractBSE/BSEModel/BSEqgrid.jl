struct BSEqgrid{T <: Union{BSEspinless, BSEspinful, BSEspinful_block}} <: AbstractBSE
	BSE::T
end
function Base.getproperty(bse::BSEqgrid, sym::Symbol)
	if sym == :BSE
		return bse.BSE
	else
		return getfield(bse.BSE, sym)
	end
end
function Base.show(io::IO, bse::BSEqgrid)
	show(io, bse.BSE)
	println("Only work for special q.")
end
"""
	struct BSEqgrid{T} <: AbstractBSE
		BSE::T
	end
	It's only a shell used for special q.
"""
function BSEqgrid(BSE::T) where {T <: Union{BSEspinless, BSEspinful, BSEspinful_block}}
	BSE.Kernal(Val(:initialize_qgrid))
	return BSEqgrid{T}(BSE)
end
function (bse::BSEqgrid{BSEspinless})(q::ReducedCoordinates)
	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)
	return bse(Htriplet, Hsinglet, q)
end
function (bse::BSEqgrid{BSEspinful})(q::ReducedCoordinates)
	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, q)
end
function (bse::BSEqgrid{BSEspinless})(Htriplet, Hsinglet, q::ReducedCoordinates)
	_BSE_preprocess_q(bse, q)
	return _BSE_Hamiltonian!(bse.BSE, Htriplet, Hsinglet)
end
function (bse::BSEqgrid{BSEspinful})(H, q::ReducedCoordinates)
	_BSE_preprocess_q!(bse, q)
	return _BSE_Hamiltonian!(bse.BSE, H)
end
function _BSE_preprocess_q!(bse::BSEqgrid, q::ReducedCoordinates)
	# We can't calculate Γ directly.
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse.BSE, q)
		bse.Kernal(q)
	else
		kgrid = bse.kgrid
		kgrid_Γ = bse.Kernal.kgrid_Γ

		nk = length(kgrid)
		kq_kindex = Vector{Int}(undef, nk)
		kΓq_kΓindex = Vector{Int}(undef, nk)
		Threads.@threads for ki in Base.OneTo(nk)
			kq = kgrid[ki] + q
			kq_kindex[ki] = findfirst(k -> all(isinteger, k - kq), kgrid)
			kΓq = kgrid_Γ[ki] + q
			kΓq_kΓindex[ki] = findfirst(k -> all(isinteger, k - kΓq), kgrid_Γ)
		end

		bse.bandkq .= bse.bandk[kq_kindex]
		bse.Kernal(kq_kindex, kΓq_kΓindex, q)
	end
	return nothing
end
function _BSE_preprocess_eleband_q!(bse::BSEqgrid, q::ReducedCoordinates)
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse.BSE, q)
	else
		kgrid = bse.kgrid
		kgrid_Γ = bse.Kernal.kgrid_Γ

		nk = length(kgrid)
		kq_kindex = Vector{Int}(undef, nk)
		kΓq_kΓindex = Vector{Int}(undef, nk)
		Threads.@threads for ki in Base.OneTo(nk)
			kq = kgrid[ki] + q
			kq_kindex[ki] = findfirst(k -> all(isinteger, k - kq), kgrid)
			kΓq = kgrid_Γ[ki] + q
			kΓq_kΓindex[ki] = findfirst(k -> all(isinteger, k - kΓq), kgrid_Γ)
		end

		bse.bandkq .= bse.bandk[kq_kindex]
	end
	return q
end
