struct ijRMap
	norb::Int
	nR::Int
	N::Int
	ijR2idx::Array{Int, 3}
	idx2ijR::Vector{Tuple{Int, Int, Int}}
end

Base.length(ijRmap::ijRMap) = ijRmap.N
Base.size(ijRmap::ijRMap) = (ijRmap.norb, ijRmap.norb, ijRmap.nR)
function Base.getindex(ijRmap::ijRMap, idx::Integer)
	return ijRmap.idx2ijR[idx]
end
function Base.getindex(ijRmap::ijRMap, ijR::Tuple{<:Integer, <:Integer, <:Integer})
	return ijRmap.ijR2idx[ijR[1], ijR[2], ijR[3]]
end
function Base.show(io::IO, ijRmap::ijRMap)
	print(io, "ijRMap(")
	print(io, "norb = ", ijRmap.norb)
	print(io, ", nR = ", ijRmap.nR)
	print(io, ")")
end

"""
	ijRMap(norb, nR)

Create a object `ijRMap`, used as the basis os exciton.
Here `norb` should be a single Integer, `nR` is also a Integer.
"""
function ijRMap(norb::Integer, nR::Integer)

	norb = Int(norb)
	nR = Int(nR)

	N = norb^2 * nR
	idx2ijR = Vector{Tuple{Int, Int, Int}}(undef, N)
	ijR2idx = Array{Int, 3}(undef, norb, norb, nR)
	n = 0
	for R in 1:nR, j in 1:norb, i in 1:norb
		n += 1
		idx2ijR[n] = (i, j, R)
		ijR2idx[i, j, R] = n
	end

	return ijRMap(norb, nR, N, ijR2idx, idx2ijR)
end
