struct ijMap
	norb::Int
	N::Int
	ij2idx::Array{Int, 2}
	idx2ij::Vector{Tuple{Int, Int}}
end

Base.length(ijmap::ijMap) = ijmap.N
Base.size(ijmap::ijMap) = (ijmap.norb, ijmap.norb)
function Base.getindex(ijmap::ijMap, idx::Integer)
	return ijmap.idx2ij[idx]
end
function Base.getindex(ijmap::ijMap, ij::Tuple{<:Integer, <:Integer})
	return ijmap.ij2idx[ij[1], ij[2]]
end
function Base.show(io::IO, ijmap::ijMap)
	print(io, "ijMap(")
	print(io, "norb = ", ijmap.norb)
	print(io, ")")
end

"""
	ijMap(norb)

Create a object `ijMap`, used as the basis os exciton.
Here `norb` should be a single Integer.
"""
function ijMap(norb::Integer)

	norb = Int(norb)

	N = norb^2
	idx2ij = Vector{Tuple{Int, Int}}(undef, N)
	ij2idx = Array{Int, 2}(undef, norb, norb)
	n = 0
	for j in 1:norb, i in 1:norb
		n += 1
		idx2ij[n] = (i, j)
		ij2idx[i, j] = n
	end

	return ijMap(norb, N, ij2idx, idx2ij)
end
