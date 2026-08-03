abstract type _Kˣ_upu end
@inline Base.getindex(upu::_Kˣ_upu, i::Int, j::Int) = upu.upu[i, j]

include("./_Kx_upu_spinaware.jl")
include("./_Kx_upu_spinblind.jl")

function _Kˣ_upu_pH(TB, kgrid)

	norb = numorb(TB)
	dorb = Matrix{eltype(TB.orb_location)}(undef, norb, norb)
	for j in 1:norb, i in 1:norb
		dorb[i, j] = TB.orb_location[i] - TB.orb_location[j]
	end

	pH = Array{ComplexF64}(undef, norb, norb, 3, length(kgrid))
	Threads.@threads for ik in eachindex(kgrid)
		k = kgrid[ik]
		pHk = view(pH,:,:,:,ik)
		TB.H(Val(:partial), pHk, TB.lattice, k, TB.orb_location)
		# later use ψk, so remove the phase factor here.
		for j in 2:norb, i in 1:(j-1)
			phase = cispi(2 * (k ⋅ dorb[i, j]))
			pHk[i, j, :] .*= phase
			pHk[j, i, :] .*= conj(phase)
		end
	end

	return pH
end
