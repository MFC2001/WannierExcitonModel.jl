
mutable struct BerryCurvature_wilsonloop
	lattice::Lattice
	center::RedKgrid
	vertex::Vector{ReducedCoordinates{Rational{Int}}}
	center_vertex_loop::Vector{Vector{Int}}
	δS::Float64
	buffer::Any
end
function (BC::BerryCurvature_wilsonloop)(u_getter::WilsonLoopWaveGetter{1})
	BC_value = similar(BC.center_vertex_loop, Float64)
	for (ki, loop) in enumerate(BC.center_vertex_loop)
		BC_value[ki] = imag(log(WilsonLoop(loop, u_getter))) / BC.δS
	end
	BC.buffer = BC_value
	return BC_value
end
function (BC::BerryCurvature_wilsonloop)(u_getter::WilsonLoopWaveGetter{N}) where {N}
	BC_value = Array{ComplexF64}(undef, N, N, length(BC.center_vertex_loop))
	for (ki, loop) in enumerate(BC.center_vertex_loop)
		BC_value[:, :, ki] .= -im * log(WilsonLoop(loop, u_getter)) / BC.δS
	end
	BC.buffer = BC_value
	return BC_value
end
function (BC::BerryCurvature_wilsonloop)(sym::Symbol, args...)
	return BC(Val(sym), args...)
end
function (BC::BerryCurvature_wilsonloop)(::Union{Val{:ChernNumber}, Val{:CN}})
	return _BerryCurvature_wilsonloop_ChernNumber(BC, BC.buffer)
end
function _BerryCurvature_wilsonloop_ChernNumber(BC::BerryCurvature_wilsonloop, BC_v::AbstractVector{<:Real})
	return sum(BC_v) * BC.δS / 2π
end
function _BerryCurvature_wilsonloop_ChernNumber(BC::BerryCurvature_wilsonloop, BC_v::AbstractArray{<:Number, 3})
	TBC = similar(BC_v, size(BC_v, 3))
	Threads.@threads for k in axes(BC_v, 3)
		TBC[k] = tr(view(BC_v, :, :, k))
	end
	return sum(TBC) * BC.δS / 2π
end
# function (BC::BerryCurvature_wilsonloop)(::Union{Val{:ChernNumber}, Val{:CN}}, BC_v::AbstractVector{<:Real})
# 	return sum(BC_v) * BC.δS / 2π
# end
# function (BC::BerryCurvature_wilsonloop)(::Union{Val{:ChernNumber}, Val{:CN}}, BC_v::AbstractArray{<:Number, 3})
# 	TBC = similar(BC_v, size(BC_v, 3))
# 	Threads.@threads for k in axes(BC_v, 3)
# 		TBC[k] = tr(view(BC_v, :, :, k))
# 	end
# 	return sum(TBC) * BC.δS / 2π
# end
function (BC::BerryCurvature_wilsonloop)(::Union{Val{:ChernNumber}, Val{:CN}}, u_getter::WilsonLoopWaveGetter{N}) where {N}
	return _BerryCurvature_wilsonloop_ChernNumber(BC, BC(u_getter))
end
"""
	BC = BerryCurvature(::Val{:WilsonLoop}, kgrid::MonkhorstPack, lattice::Lattice, kz::Real = 0)
	kgrid should only have one component is 1, 
	and kz represents the reduced coordinates at the same direction as the direction where the component of kgrid is 1. 

	BC is a instance of BerryCurvature_wilsonloop, 
	it have vertex field, which is used to calculate the periodic part of bloch wave function.
	Then create an instance u_getter of WilsonLoopWaveGetter{N}, BC(ugetter) will return Berry Curvature.
"""
function BerryCurvature(::Union{Val{:WilsonLoop}, Val{:WL}}, kgrid::MonkhorstPack, lattice::Lattice, kz::Real = 0)

	center, vertex, center_vertex_loop = _BerryCurvature_center_vertex(kgrid.kgrid_size, kz)
	δS = _BerryCurvature_microsurface(kgrid.kgrid_size, lattice)

	# fold all kpoints into the first Brillouin zone
	center = fold2FBZ(reciprocal(lattice), center)
	kshift = map(nk -> isodd(nk) ? 0 : 1 // 2, kgrid.kgrid_size)
	center = RedKgrid(center, kgrid.kgrid_size, kshift)

	return BerryCurvature_wilsonloop(lattice, center, vertex, center_vertex_loop, δS, nothing)
end
function _BerryCurvature_center_vertex(kgrid_size, kz::T) where {T <: Real}
	count(kgrid_size .== 1) == 1 || error("Only 2D kgrid are supported.")

	d1 = findfirst(kgrid_size .== 1)
	if d1 == 1
		Nkx = kgrid_size[2]
		Nky = kgrid_size[3]
		Ixy = [2, 3]
		Iz = 1
	elseif d1 == 2
		Nkx = kgrid_size[3]
		Nky = kgrid_size[1]
		Ixy = [3, 1]
		Iz = 2
	else
		Nkx = kgrid_size[1]
		Nky = kgrid_size[2]
		Ixy = [1, 2]
		Iz = 3
	end

	Δkx = 1 // Nkx
	Δky = 1 // Nky

	Tprom = promote_type(T, Rational{Int})

	vertex_x = map(1:Nkx+1) do nk
		(nk - 1) * Δkx - 1 // 2
	end
	vertex_y = map(1:Nky+1) do nk
		(nk - 1) * Δky - 1 // 2
	end
	vertex = vec([(x, y) for x in vertex_x, y in vertex_y])
	vertex = map(vertex) do v
		Tv = Tprom[0, 0, 0]
		Tv[Ixy] .= v
		Tv[Iz] = kz
		return ReducedCoordinates(Tv)
	end

	center_x = map(1:Nkx) do nk
		(nk - 1 // 2) * Δkx - 1 // 2
	end
	center_y = map(1:Nky) do nk
		(nk - 1 // 2) * Δky - 1 // 2
	end
	center = vec([(x, y) for x in center_x, y in center_y])
	center = map(center) do c
		Tc = Tprom[0, 0, 0]
		Tc[Ixy] .= c
		Tc[Iz] = kz
		return ReducedCoordinates(Tc)
	end

	center_vertex_loop = Matrix{Vector{Tuple{Int, Int}}}(undef, Nkx, Nky)
	for xi in 1:Nkx, yi in 1:Nky
		center_vertex_loop[xi, yi] = [
			(xi, yi),
			(xi + 1, yi),
			(xi + 1, yi + 1),
			(xi, yi + 1),
			(xi, yi),
		]
	end
	center_vertex_loop = reshape(center_vertex_loop, :)
	center_vertex_loop = map(loop -> map(v -> v[1] + (v[2] - 1) * (Nkx + 1), loop), center_vertex_loop)

	return center, vertex, center_vertex_loop
end
function _BerryCurvature_microsurface(kgrid_size, lattice)
	count(kgrid_size .== 1) == 1 || error("Only 2D kgrids are supported.")
	𝐚, 𝐛, 𝐜 = basisvectors(reciprocal(lattice))
	if kgrid_size[1] == 1
		S = abs((𝐛×𝐜)[3])
	elseif kgrid_size[2] == 1
		S = abs((𝐚×𝐜)[3])
	else
		S = abs((𝐚×𝐛)[3])
	end
	return S / prod(kgrid_size)
end
