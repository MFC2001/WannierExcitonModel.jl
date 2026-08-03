struct _Kˣ_upu_spinblind <: _Kˣ_upu
	upu::Matrix{ComplexF64} # [1:3, ivck]
	CoulombScaledivnk::Float64
	rlattice::ReciprocalLattice{Float64}
end
function _Kˣ_upu_spinblind(TB::AbstractTightBindModel, kgrid, bandk, vckmap)
	np = count(TB.period)
	if np ≠ 3
		upu = Matrix{ComplexF64}(undef, 0, 0)
		return _Kˣ_upu_spinblind(upu, 0.0, reciprocal(TB.lattice))
	end

	pH = _Kˣ_upu_pH(TB, kgrid)

	upu = Matrix{ComplexF64}(undef, 3, length(vckmap))
	Threads.@threads for i in Base.OneTo(length(vckmap))
		(v, c, k) = vckmap[i]
		ψv = view(bandk[k].vectors, :, v)
		ψc = view(bandk[k].vectors, :, c)
		ΔE = bandk[k].values[c] - bandk[k].values[v]
		for α in 1:3
			pHk = view(pH,:,:,α,k)
			upu[α, i] = (ψv ⋅ (pHk * ψc)) / ΔE
		end
	end

	Ω = abs(det(parent(TB.lattice)))
	CoulombScaledivnk = CoulombScale * 4π / Ω / length(kgrid)

	rlattice = reciprocal(TB.lattice)

	return _Kˣ_upu_spinblind(upu, CoulombScaledivnk, rlattice)
end
