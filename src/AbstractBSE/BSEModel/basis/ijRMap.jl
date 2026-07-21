struct ijRMap
	nw_h::Int
	nw_e::Int
	nR::Int
	N::Int
	ijR2idx::Array{Int, 3}
	idx2ijR::Vector{Tuple{Int, Int, Int}}
end
Base.iterate(ijRmap::ijRMap, state = 1) = state > length(ijRmap) ? nothing : (ijRmap[state], state + 1)
Base.eltype(::ijRMap) = Tuple{Int, Int, Int}
Base.length(ijRmap::ijRMap) = ijRmap.N
Base.size(ijRmap::ijRMap) = (ijRmap.nw_h, ijRmap.nw_e, ijRmap.nR)
Base.getindex(ijRmap::ijRMap, idx::Integer) = ijRmap.idx2ijR[idx]
Base.getindex(ijRmap::ijRMap, ijR::Tuple{<:Integer, <:Integer, <:Integer}) =
	ijRmap.ijR2idx[ijR[1], ijR[2], ijR[3]]
Base.firstindex(::ijRMap) = 1
Base.lastindex(ijRmap::ijRMap) = ijRmap.N
function Base.show(io::IO, ijRmap::ijRMap)
	print(io, "ijRMap(")
	print(io, "nw_h = ", ijRmap.nw_h)
	print(io, ", nw_e = ", ijRmap.nw_e)
	print(io, ", nR = ", ijRmap.nR)
	print(io, ")")
end

"""
	ijRMap(nw_h::Integer, nw_e::Integer, nR::Integer)

Create a object `ijRMap`, used as the basis of exciton.

	ijRMap(nw::Integer, nR::Integer)

this means `nw_h` equal to `nw_e`.
"""
ijRMap(nw::Integer, nR::Integer) = ijRMap(nw, nw, nR)
function ijRMap(nw_h::Integer, nw_e::Integer, nR::Integer)

	nw_h = Int(nw_h)
	nw_e = Int(nw_e)
	nR = Int(nR)

	N = nw_h * nw_e * nR
	idx2ijR = Vector{Tuple{Int, Int, Int}}(undef, N)
	ijR2idx = Array{Int, 3}(undef, nw_h, nw_e, nR)
	n = 0
	for R in 1:nR, j in 1:nw_e, i in 1:nw_h
		n += 1
		idx2ijR[n] = (i, j, R)
		ijR2idx[i, j, R] = n
	end

	return ijRMap(nw_h, nw_e, nR, N, ijR2idx, idx2ijR)
end
