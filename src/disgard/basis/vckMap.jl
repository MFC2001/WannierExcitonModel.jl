"""
    vckMap(v, c, nk) -> vckMap

Construct a `vckMap` object representing the basis of exciton transitions.

**Arguments**:
- `v`: Valence band index (or vector of indices).
- `c`: Conduction band index (or vector of indices).
- `nk`: Number of k‑points (integer).

The resulting object has size `(length(v), length(c), nk)` and supports two‑way indexing:

```julia
julia> vckmap = vckMap([1, 2], [3, 4], 9)
vckMap(v = [1, 2], c = [3, 4], nk = 9)

julia> length(vckmap) == 2 * 2 * 9
true

julia> vckmap[1]         # linear index → tuple (v, c, k)
(1, 3, 1)

julia> vckmap[1, 3, 1]   # tuple → linear index
1
```
"""
struct vckMap
    nv::Int
    nc::Int
    nk::Int
    N::Int
    v2idx::Dict{Int,Int}
    c2idx::Dict{Int,Int}
    idx2v::Vector{Int}
    idx2c::Vector{Int}
    vck2idx::Array{Int,3}
    idx2vck::Vector{Tuple{Int,Int,Int}}
end
function Base.getproperty(vckmap::vckMap, name::Symbol)
    if name === :v
        return getfield(vckmap, :idx2v)
    elseif name === :c
        return getfield(vckmap, :idx2c)
    else
        return getfield(vckmap, name)
    end
end
Base.iterate(vckmap::vckMap, state=1) = state > length(vckmap) ? nothing : (vckmap[state], state + 1)
Base.eltype(::vckMap) = Tuple{Int,Int,Int}
Base.length(vckmap::vckMap) = vckmap.N
Base.size(vckmap::vckMap) = (vckmap.nv, vckmap.nc, vckmap.nk)
Base.getindex(vckmap::vckMap, idx::Integer) = vckmap.idx2vck[idx]
Base.getindex(vckmap::vckMap, vck::Tuple{<:Integer,<:Integer,<:Integer}) =
    vckmap.vck2idx[vckmap.v2idx[vck[1]], vckmap.c2idx[vck[2]], vck[3]]
Base.firstindex(::vckMap) = 1
Base.lastindex(vckmap::vckMap) = vckmap.N
function Base.show(io::IO, vckmap::vckMap)
    print(io, "vckMap(")
    print(io, "v = ", vckmap.idx2v)
    print(io, ", c = ", vckmap.idx2c)
    print(io, ", nk = ", vckmap.nk)
    print(io, ")")
end
function Base.reshape(A::AbstractVector, vckmap::vckMap)
    @assert length(A) == length(vckmap) "Vector length mismatch"
    return reshape(A, vckmap.nv, vckmap.nc, vckmap.nk)
end
function Base.reshape(A::AbstractArray{T,3}, vckmap::vckMap) where {T}
    @assert size(A) == size(vckmap) "Array size mismatch"
    return reshape(A, :)
end
"""
	vckMap(v, c, nk)

Create a object `vckMap`, used as the basis os exciton.
Here `v` and `c` can be a vector of Integer or a single Integer, `nk` is a Integer.
"""
function vckMap(v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, nk::Integer)

    idx2v = map(Int, collect(v))
    idx2c = map(Int, collect(c))

    nv = length(idx2v)
    nc = length(idx2c)
    nk = Int(nk)

    v2idx = Dict{Int,Int}()
    for i in eachindex(idx2v)
        v2idx[idx2v[i]] = i
    end
    c2idx = Dict{Int,Int}()
    for i in eachindex(idx2c)
        c2idx[idx2c[i]] = i
    end

    N = nv * nc * nk
    idx2vck = Vector{Tuple{Int,Int,Int}}(undef, N)
    vck2idx = Array{Int,3}(undef, nv, nc, nk)
    n = 0
    for k in 1:nk, c in 1:nc, v in 1:nv
        n += 1
        idx2vck[n] = (idx2v[v], idx2c[c], k)
        vck2idx[v, c, k] = n
    end

    return vckMap(nv, nc, nk, N, v2idx, c2idx, idx2v, idx2c, vck2idx, idx2vck)
end
vckMap(v::AbstractVector{<:Integer}, c::Integer, nk::Integer) = vckMap(v, [c], nk)
vckMap(v::Integer, c::AbstractVector{<:Integer}, nk::Integer) = vckMap([v], c, nk)
vckMap(v::Integer, c::Integer, nk::Integer) = vckMap([v], [c], nk)




# struct vckVector{T <: Number} <: AbstractVector{T}
# 	data::Array{T, 3}
# 	vckmap::vckMap
# end
# struct vckMatrix{T <: Number} <: AbstractMatrix{T}
# 	data::Array{T, 6}
# 	vckmap::vckMap
# end


# function vckVector(
# 	v_vals::AbstractVector{<:Integer},
# 	c_vals::AbstractVector{<:Integer},
# 	nk::Integer,
# 	data = Array{T, 3}(undef, length(v_vals), length(c_vals), nk),
# )

# end

# Base.size(H::BSEHamiltonian) = (H.nv, H.nc, H.nk, H.nv, H.nc, H.nk)
# Base.getindex(H::HermitianStorage, i::Int, j::Int) = begin
# 	if H.part == :U
# 		i ≤ j ? H.storage[triu_index(i, j, H.n)] : conj(H.storage[triu_index(j, i, H.n)])
# 	else
# 		i ≥ j ? H.storage[tril_index(i, j, H.n)] : conj(H.storage[tril_index(j, i, H.n)])
# 	end
# end
