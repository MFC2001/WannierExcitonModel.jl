
function ChernNumber(sym::Symbol, paras...; kwargs...)
	return ChernNumber(Val(sym), paras...; kwargs...)
end
"""
	ChernNumber(BC, kgrid::AbstractBrillouinZone, lattice::Lattice)
	BC is a Vector{<:Real} or Array{<:Number, 3}, Vector{<:Real} means abelian, Array{<:Number, 3} means non-abelian.
	when non-abelian, size(BC) = (nband, nband, nk), CN = tr(BC).
"""
function ChernNumber(BC::AbstractVector{<:Real}, kgrid::AbstractBrillouinZone, lattice::Lattice)
	𝐚, 𝐛, 𝐜 = basisvectors(reciprocal(lattice))
	if kgrid.kgrid_size[1] == 1
		S = abs((𝐛×𝐜)[3])
	elseif kgrid.kgrid_size[2] == 1
		S = abs((𝐚×𝐜)[3])
	elseif kgrid.kgrid_size[3] == 1
		S = abs((𝐚×𝐛)[3])
	else
		error("Only calculate Chern Number for 2D.")
	end

	return sum(BC) * S / length(kgrid) / 2π
end
function ChernNumber(BC::AbstractArray{<:Number, 3}, kgrid::AbstractBrillouinZone, lattice::Lattice)
	𝐚, 𝐛, 𝐜 = basisvectors(reciprocal(lattice))
	if kgrid.kgrid_size[1] == 1
		S = abs((𝐛×𝐜)[3])
	elseif kgrid.kgrid_size[2] == 1
		S = abs((𝐚×𝐜)[3])
	elseif kgrid.kgrid_size[3] == 1
		S = abs((𝐚×𝐛)[3])
	else
		error("Only calculate Chern Number for 2D.")
	end

	TBC = similar(BC, size(BC, 3))
	Threads.@threads for k in axes(BC, 3)
		TBC[k] = tr(BC[:, :, k])
	end
	return sum(TBC) * S / length(kgrid) / 2π
end