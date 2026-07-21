struct _Kˣ_upu_spinaware <: _Kˣ_upu
	upu_uu::Matrix{ComplexF64}
	upu_dd::Matrix{ComplexF64}
	CoulombScaledivnk::Float64
	rlattice::ReciprocalLattice{Float64}
end
function _Kˣ_upu_spinaware(TB_up, TB_dn, kgrid, bandk_up, bandk_dn, vckmap_uu, vckmap_dd)
	np = count(TB_up.period)
	if np ≠ 3
		upu_uu = Matrix{ComplexF64}(undef, 0, 0)
		upu_dd = Matrix{ComplexF64}(undef, 0, 0)
		return _Kˣ_upu_spinaware(upu_uu, upu_dd, 0.0, reciprocal(TB_up.lattice))
	end

	pH_up = _Kˣ_upu_pH(TB_up, kgrid)
	pH_dn = _Kˣ_upu_pH(TB_dn, kgrid)

	n_uu = length(vckmap_uu)
	upu_uu = Matrix{ComplexF64}(undef, 3, n_uu)
	if n_uu > 0
		Threads.@threads for i in Base.OneTo(n_uu)
			(v, c, k) = vckmap_uu[i]
			ψv = view(bandk_up[k].vectors, :, v)
			ψc = view(bandk_up[k].vectors, :, c)
			ΔE = bandk_up[k].values[c] - bandk_up[k].values[v]
			for α in 1:3
				pHk = view(pH_up,:,:,α,k)
				upu_uu[α, i] = (ψv ⋅ (pHk * ψc)) / ΔE
			end
		end
	end

	n_dd = length(vckmap_dd)
	upu_dd = Matrix{ComplexF64}(undef, 3, n_dd)
	if n_dd > 0
		Threads.@threads for i in Base.OneTo(n_dd)
			(v, c, k) = vckmap_dd[i]
			ψv = view(bandk_dn[k].vectors, :, v)
			ψc = view(bandk_dn[k].vectors, :, c)
			ΔE = bandk_dn[k].values[c] - bandk_dn[k].values[v]
			for α in 1:3
				pHk = view(pH_dn,:,:,α,k)
				upu_dd[α, i] = (ψv ⋅ (pHk * ψc)) / ΔE
			end
		end
	end

	Ω = abs(det(parent(TB_up.lattice)))
	CoulombScaledivnk = CoulombScale * 4π / Ω / length(kgrid)

	rlattice = reciprocal(TB_up.lattice)

	return _Kˣ_upu_spinaware(upu_uu, upu_dd, CoulombScaledivnk, rlattice)
end
