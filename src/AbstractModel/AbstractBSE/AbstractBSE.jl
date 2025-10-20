
export BSE, BSEwannier

abstract type AbstractBSE <: AbstractModel end

include("./AbstractKernalInterAction/AbstractKernalInterAction.jl")
include("./others/others.jl")

include("./BSEModel/BSEModel.jl")


"""
	BSE(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction, sym::Symbol;
		kgrid::Union{MonkhorstPack, RedKgrid} = MonkhorstPack([1, 1, 1]), 
		v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, 
		scissor::Real = 0, 
		isqgrid::Bool = false) -> AbstractBSE

Create a BSE model, the type depends on `sym`.
- `TB::AbstractTightBindModel`: contains crystal structure and electronic band structure;
- `Kernal::AbstractKernalInterAction`: describes how to calculate `Kᵈ` and `Kˣ`;
- `sym::Symbol`: can be set as `:SU2`, `:general`, `:cluster_SU2`, `:cluster_general`;
- `kgrid::Union{MonkhorstPack, RedKgrid}`: defines the electronic kgrid used to construct the basis of exciton;
- `v::AbstractVector{<:Integer}`: the valence band index used to construct the basis of exciton;
- `c::AbstractVector{<:Integer}`: the conduction band index used to construct the basis of exciton;
- `scissor::Real`: scissor operator used to provide electronic band gap correction;
- `isqgrid::Bool`: when all the excitonic momentum needed can be regarded as the momentum difference in the kgrid.
"""
function BSE(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction, sym::Symbol;
	kgrid::Union{MonkhorstPack, RedKgrid} = MonkhorstPack([1, 1, 1]), v, c, scissor::Real = 0, isqgrid::Bool = false)
	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	end
	return BSE(Val(sym), TB, Kernal; kgrid, v, c, scissor, isqgrid)
end
function BSE(::Val{:SU2}, TB, Kernal; kwards...)
	return BSESU2(TB, Kernal; kwards...)
end
function BSE(::Val{:general}, TB, Kernal; kwards...)
	return BSEgeneral(TB, Kernal; kwards...)
end
function BSE(::Val{:general_block}, TB, Kernal; kwards...)
	error("To be continued.")
	return BSEgeneral_block(TB, Kernal; kwards...)
end
function BSE(::Val{:cluster_SU2}, TB, Kernal; kwards...)
	return BSEcluster_SU2(TB, Kernal; kwards...)
end
function BSE(::Val{:cluster_general}, TB, Kernal; kwards...)
	return BSEcluster_general(TB, Kernal; kwards...)
end
function BSE(::Val{:cluster_general_block}, TB, Kernal; kwards...)
	error("To be continued.")
	return BSEcluster_general_block(TB, Kernal; kwards...)
end
"""
	BSEwannier(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction, sym::Symbol;
		qgrid::Union{MonkhorstPack, RedKgrid}, v, c, scissor::Real = 0)

Create a BSE used to calculate excitonic maximally localized wannier function.
This function is only a shortcut, it will help you set:
- `kgrid`: is twice as much as `qgrid`;
- `isqgrid`: set as `true`.
Then you can run:
```julia
BSE_wannier(qgrid, bse; kwargs...)
```
"""
function BSEwannier(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction, sym::Symbol;
	qgrid::Union{MonkhorstPack, RedKgrid}, v, c, scissor::Real = 0)

	p = count(TB.period)
	if p == 3
		kgrid = qgrid.kgrid_size .* 2
		kshift = [1 // 2, 1 // 2, 1 // 2]
	elseif p == 2
		kgrid = collect(qgrid.kgrid_size .* 2)
		kshift = [1 // 2, 1 // 2, 1 // 2]
		if !TB.period[1]
			kgrid[1] = 1
			kshift[1] = 0
		elseif !TB.period[2]
			kgrid[2] = 1
			kshift[2] = 0
		elseif !TB.period[3]
			kgrid[3] = 1
			kshift[3] = 0
		end
	elseif p == 1
		if TB.period[1]
			kgrid = [qgrid.kgrid_size[1] * 2, 1, 1]
			kshift = [1 // 2, 0, 0]
		elseif TB.period[2]
			kgrid = [1, qgrid.kgrid_size[1] * 2, 1]
			kshift = [0, 1 // 2, 0]
		elseif TB.period[3]
			kgrid = [1, 1, qgrid.kgrid_size[1] * 2]
			kshift = [0, 0, 1 // 2]
		end
	elseif p == 0
		kgrid = [1, 1, 1]
		kshift = [0, 0, 0]
	end
	#basis is twice as dense
	kgrid = MonkhorstPack(kgrid; kshift)

	return BSE(TB, Kernal, sym; kgrid, v, c, scissor, isqgrid = true)
end


# function Base.setproperty!(bse::BSE, prop::Symbol, val)
# 	if prop == :kgrid
# 		if val isa MonkhorstPack
# 			val = RedKgrid(val)
# 		elseif val isa RedKgrid
# 		else
# 			error("Wrong kgrid to BSE.")
# 		end
# 		setfield!(bse, prop, val)
# 		vckmap = vckMap(bse.vckmap.idx2v, bse.vckmap.idx2c, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 		bandk = BAND(val, bse.TB; vector = true)
# 		setfield!(bse, :bandk, bandk)
# 	elseif prop == :v
# 		vckmap = vckMap(prop, bse.vckmap.idx2c, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 	elseif prop == :c
# 		vckmap = vckMap(bse.vckmap.idx2v, prop, length(bse.kgrid))
# 		setfield!(bse, :vckmap, vckmap)
# 	elseif prop == :vckmap
# 		if prop.nk ≠ length(bse.kgrid)
# 			error("Please check the vckmap.")
# 		end
# 		setfield!(bse, :vckmap, prop)
# 	elseif prop == :bandk
# 		error("Please don't change BSE.vckmap or BSE.bandk directly.")
# 	else
# 		setfield!(bse, prop, val)
# 	end
# end




