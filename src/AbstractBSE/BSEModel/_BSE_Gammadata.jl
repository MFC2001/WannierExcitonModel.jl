function _pH_TB_kgrid(TB, kgrid)

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

struct _BSE_Γdata
	upHu::Matrix{ComplexF64} # [1:3, ivck], eV ⋅ Å
	upu::Matrix{ComplexF64} # [1:3, ivck], Å
	CoulombScaledivnk::Float64
	rlattice::ReciprocalLattice{Float64}
end
function _BSE_Γdata(TB::AbstractTightBindModel, kgrid, bandk, vckmap)

	Ω = abs(det(parent(TB.lattice)))
	CoulombScaledivnk = CoulombScale * 4π / Ω / length(kgrid)

	rlattice = reciprocal(TB.lattice)

	Nvck = length(vckmap)
	if Nvck == 0
		upHu = Matrix{ComplexF64}(undef, 3, 0)
		upu = Matrix{ComplexF64}(undef, 3, 0)
		return _BSE_Γdata(upHu, upu, CoulombScaledivnk, rlattice)
	end

	pH = _pH_TB_kgrid(TB, kgrid)

	upHu = Matrix{ComplexF64}(undef, 3, Nvck)
	Threads.@threads for i in Base.OneTo(Nvck)
		(v, c, k) = vckmap[i]
		ψv = view(bandk[k].vectors, :, v)
		ψc = view(bandk[k].vectors, :, c)
		for α in 1:3
			pHk = view(pH,:,:,α,k)
			upHu[α, i] = ψv ⋅ (pHk * ψc) # eV ⋅ Å
		end
	end

	if count(TB.period) == 3
		upu = Matrix{ComplexF64}(undef, 3, Nvck) # Å
		Threads.@threads for i in Base.OneTo(Nvck)
			(v, c, k) = vckmap[i]
			ΔE = bandk[k].values[c] - bandk[k].values[v] # eV
			upu[:, i] .= upHu[:, i] ./ ΔE
		end
	else
		upu = Matrix{ComplexF64}(undef, 3, 0)
	end

	return _BSE_Γdata(upHu, upu, CoulombScaledivnk, rlattice)
end
