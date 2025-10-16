
export MirrorCorrection

include("./MirrorCorrection_2D.jl")
include("./MirrorCorrection_3D.jl")

"""
	MirrorCorrection(U::HR, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates}, rcut::Real; kwargs...) -> HR

Try to apply a mirror correction to `U`.
This mirror correction is achieved by using Gauss potential.

- `rcut::Real`: should be less than a real radius value, out of which the U approximate to 1/r, default is judged by `U` and `TB`;
kwargs:
- `period`: the periodicity;
- `kgrid::MonkhorstPack`: should be the kgrid used to calculate `U`, default is judged by `U`;
- `ϵ::Real`: dielectric constant;
- `αrcut::Real`: α × rcut, default is 4.5;
- `δ:Real`: the cutoff value of reciprocal potential.
"""
function MirrorCorrection(U::HR, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates}, rcut::Union{Nothing, Real} = nothing;
	period = Bool[1, 1, 1], kgrid = nothing, αrcut = 4.5, δ = 1e-8, ϵ = 1, head = nothing)

	if isnothing(kgrid)
		kgrid_max = map(maximum, eachcol(U.path[:, 1:3]))
		kgrid_min = map(minimum, eachcol(U.path[:, 1:3]))
		kgrid = MonkhorstPack((kgrid_max - kgrid_min) .+ 1)
		kgrid = RedKgrid(kgrid)
	elseif kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	else
		error("Wrong kgrid.")
	end

	p = count(period)
	if p == 3
		if isnothing(rcut)
			a₁ = lattice[:, 1]
			a₂ = lattice[:, 2]
			a₃ = lattice[:, 3]
			Ω = abs((a₁ × a₂) ⋅ a₃)
			h₁ = Ω / norm(a₂ × a₃)
			h₂ = Ω / norm(a₃ × a₁)
			h₃ = Ω / norm(a₁ × a₂)
			half_size = floor.(kgrid.kgrid_size / 2)
			rcut = min(h₁ * half_size[1], h₂ * half_size[2], h₃ * half_size[3])
		end
		return MirrorCorrection_3D(U, kgrid, lattice, orblocat, rcut; αrcut, δ, ϵ, head)
	elseif p == 2
		if isnothing(rcut)
			half_size = floor.(kgrid.kgrid_size / 2)
			if !period[1]
				a₂ = lattice[:, 2]
				a₃ = lattice[:, 3]
				S = abs((a₂ × a₃))
				h₂ = S / norm(a₃)
				h₃ = S / norm(a₂)
				rcut = min(h₂ * half_size[2], h₃ * half_size[3])
			elseif !period[2]
				a₃ = lattice[:, 3]
				a₁ = lattice[:, 1]
				S = abs((a₃ × a₁))
				h₃ = S / norm(a₁)
				h₁ = S / norm(a₃)
				rcut = min(h₁ * half_size[1], h₃ * half_size[3])
			else
				a₁ = lattice[:, 1]
				a₂ = lattice[:, 2]
				S = abs((a₁ × a₂))
				h₁ = S / norm(a₂)
				h₂ = S / norm(a₁)
				rcut = min(h₁ * half_size[1], h₂ * half_size[2])
			end
		end
		return MirrorCorrection_2D(U, kgrid, lattice, orblocat, rcut; αrcut, δ, ϵ, head, period = period)
	elseif p == 1
		#TODO
		error("To be continued.")
	else
		error("To be continued.")
	end
end
