struct vcMap
	nv::Int
	nc::Int
	N::Int
	v2idx::Dict{Int, Int}
	c2idx::Dict{Int, Int}
	idx2v::Vector{Int}
	idx2c::Vector{Int}
	vc2idx::Array{Int, 2}
	idx2vc::Vector{Tuple{Int, Int}}
end

Base.length(vcmap::vcMap) = vcmap.N
Base.size(vcmap::vcMap) = (vcmap.nv, vcmap.nc)
function Base.getindex(vcmap::vcMap, idx::Integer)
	return vcmap.idx2vc[idx]
end
function Base.getindex(vcmap::vcMap, vc::Tuple{<:Integer, <:Integer})
	return vcmap.vc2idx[vcmap.v2idx[vc[1]], vcmap.c2idx[vc[2]]]
end
function Base.show(io::IO, vcmap::vcMap)
	print(io, "vcMap(")
	print(io, "v = ", vcmap.idx2v)
	print(io, ", c = ", vcmap.idx2c)
	print(io, ")")
end

"""
	vcMap(v, c)

Create a object `vcMap`, used as the basis os exciton.
Here `v` and `c` can be a vector of Integer or a single Integer.
"""
function vcMap(v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer})

	idx2v = map(Int, v)
	idx2c = map(Int, c)

	nv = length(idx2v)
	nc = length(idx2c)

	v2idx = Dict{Int, Int}()
	for i in eachindex(idx2v)
		v2idx[idx2v[i]] = i
	end
	c2idx = Dict{Int, Int}()
	for i in eachindex(idx2c)
		c2idx[idx2c[i]] = i
	end

	N = nv * nc
	idx2vc = Vector{Tuple{Int, Int}}(undef, N)
	vc2idx = Array{Int, 2}(undef, nv, nc)
	n = 0
	for c in 1:nc, v in 1:nv
		n += 1
		idx2vc[n] = (idx2v[v], idx2c[c])
		vc2idx[v, c] = n
	end

	return vcMap(nv, nc, N, v2idx, c2idx, idx2v, idx2c, vck2idx, idx2vck)
end
vcMap(v::AbstractVector{<:Integer}, c::Integer) = vcMap(v, [c])
vcMap(v::Integer, c::AbstractVector{<:Integer}) = vcMap([v], c)
vcMap(v::Integer, c::Integer) = vcMap([v], [c])
