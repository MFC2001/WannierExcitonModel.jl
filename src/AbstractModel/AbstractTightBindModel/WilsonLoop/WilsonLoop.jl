export WilsonLoop
function WilsonLoop(TB::AbstractTightBindModel, bandindex::Integer = 1; Nkx = 51, Nky = 51) where {T}

	start = -floor.(Int, ([Nkx, Nky] .- 1) .// 2)
	stop = ceil.(Int, ([Nkx, Nky] .- 1) .// 2)

	#W(ky)
	ky = start[2]-1:stop[2]
	θ_ky = Vector{Float64}(undef, length(ky))
	Threads.@threads for kyi in eachindex(ky)
		aimkdirects = [Vec3([kx // Nkx, ky[kyi] // Nky, 0]) for kx in start[1]:stop[1]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nkx-1
			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
			W *= uu
		end
		uu = sum(i -> cis(2π * ([1, 0, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
		W = uu * W

		θ_ky[kyi] = angle(W)
	end
	θ_ky = smooth_AngleSingular(θ_ky, π)
	wcc_ky = θ_ky ./ 2π

	#W(kx)
	kx = start[1]-1:stop[1]
	θ_kx = Vector{Float64}(undef, length(kx))
	Threads.@threads for kxi in eachindex(kx)
		aimkdirects = [Vec3([kx[kxi] // Nkx, ky // Nky, 0]) for ky in start[2]:stop[2]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nky-1
			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
			W *= uu
		end
		uu = sum(i -> cis(2π * ([0, 1, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
		W = uu * W

		θ_kx[kxi] = angle(W)
	end
	θ_kx = smooth_AngleSingular(θ_kx, π)
	wcc_kx = θ_kx ./ 2π

	return wcc_kx, wcc_ky
end

function WilsonLoop(TB::AbstractTightBindModel, bandindex::AbstractVector{<:Integer} = [1, 2]; Nkx = 51, Nky = 51) where {T}

	start = -floor.(Int, ([Nkx, Nky] .- 1) .// 2)
	stop = ceil.(Int, ([Nkx, Nky] .- 1) .// 2)

	nband = length(bandindex)

	#W(ky)
	ky = start[2]-1:stop[2]
	θ_ky = Matrix{Float64}(undef, nband, length(ky))
	Threads.@threads for kyi in eachindex(ky)
		aimkdirects = [Vec3([kx // Nkx, ky[kyi] // Nky, 0]) for kx in start[1]:stop[1]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nkx-1
			uu = [uband[i+1].vectors[:, m] ⋅ uband[i].vectors[:, n] for m in bandindex, n in bandindex]
			F = svd(uu)
			uu = F.U * F.Vt
			W = uu * W
		end
		uu = [sum(i -> cis(2π * ([1, 0, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, m]) * uband[end].vectors[i, n], eachindex(TB.orb_location)) for m in bandindex, n in bandindex]
		W = uu * W

		θ_ky[:, kyi] = eigvals!(Hermitian(W))
	end
	wcc_ky = θ_ky ./ 2π


	#W(kx)
	kx = start[1]-1:stop[1]
	θ_kx = Matrix{Float64}(undef, nband, length(kx))
	Threads.@threads for kxi in eachindex(kx)
		aimkdirects = [Vec3([kx[kxi] // Nkx, ky // Nky, 0]) for ky in start[2]:stop[2]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nky-1
			uu = [uband[i+1].vectors[:, m] ⋅ uband[i].vectors[:, n] for m in bandindex, n in bandindex]
			F = svd(uu)
			uu = F.U * F.Vt
			W = uu * W
		end
		uu = [sum(i -> cis(2π * ([0, 1, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, m]) * uband[end].vectors[i, n], eachindex(TB.orb_location)) for m in bandindex, n in bandindex]
		W = uu * W

		θ_kx[kxi] = eigvals!(Hermitian(W))
	end
	wcc_kx = θ_kx ./ 2π

	return wcc_kx, wcc_ky
end


function WilsonLoop(::Val{:x}, kgrid::MonkhorstPack, TB::AbstractTightBindModel, bandindex::Integer = 1)

	start = -floor.(Int, ([Nkx, Nky] .- 1) .// 2)
	stop = ceil.(Int, ([Nkx, Nky] .- 1) .// 2)

	#W(kx)
	kx = start[1]-1:stop[1]
	θ_kx = Vector{Float64}(undef, length(kx))
	Threads.@threads for kxi in eachindex(kx)
		aimkdirects = [Vec3([kx[kxi] // Nkx, ky // Nky, 0]) for ky in start[2]:stop[2]]

		uband = BAND(aimkdirects, TB, TB.orb_location; vector = true)

		W = 1
		for i in 1:Nky-1
			uu = uband[i+1].vectors[:, bandindex] ⋅ uband[i].vectors[:, bandindex]
			W *= uu
		end
		uu = sum(i -> cis(2π * ([0, 1, 0] ⋅ TB.orb_location[i])) * conj(uband[1].vectors[i, bandindex]) * uband[end].vectors[i, bandindex], eachindex(TB.orb_location))
		W = uu * W

		θ_kx[kxi] = angle(W)
	end
	θ_kx = smooth_AngleSingular(θ_kx, π)
	wcc_kx = θ_kx ./ 2π

	return wcc_kx, wcc_ky
end


function _WilsonLoop_tightbinding(kpoints, hrh, orb_location, bandindex::Integer)
	u1 = eigen!(hrh(kpoints[1], orb_location)).vectors[:, bandindex]
	u0 = u1
	W = 1
	for i in 2:length(kpoints)
		u2 = eigen!(hrh(kpoints[i], orb_location)).vectors[:, bandindex]
		uu = u2 ⋅ u1
		W *= uu
		u1 = u2
	end
	G = round.(Int, kpoints[end] - kpoints[1])
	u0 .*= cis.(2π .* (G .⋅ orb_location))
	W *= u0 ⋅ u1
	return W
end
function _WilsonLoop_tightbinding(kpoints, hrh, orb_location, bandindex::AbstractVector{<:Integer})
	u1 = eigen!(hrh(kpoints[1], orb_location)).vectors[:, bandindex]
	u0 = u1
	nband = length(bandindex)
	W = 1
	for i in 2:length(kpoints)
		u2 = eigen!(hrh(kpoints[i], orb_location)).vectors[:, bandindex]
		uu = [u2[:, m] ⋅ u1[:, n] for m in Base.OneTo(nband), n in Base.OneTo(nband)]
		F = svd(uu)
		uu = F.U * F.Vt
		W = uu * W
		u1 = u2
	end
	G = round.(Int, kpoints[end] - kpoints[1])
	u0 .*= vec(cis.(2π .* (G .⋅ orb_location)))
	uu = [u0[:, m] ⋅ u1[:, n] for m in Base.OneTo(nband), n in Base.OneTo(nband)]
	F = svd(uu)
	uu = F.U * F.Vt
	W = uu * W
	return W
end
