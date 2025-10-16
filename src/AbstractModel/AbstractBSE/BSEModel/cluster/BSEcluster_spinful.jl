struct BSEcluster_spinful{
	TBT <: AbstractTightBindModel,
	KT <: AbstractKernalInterAction,
} <: AbstractBSE
	TB::TBT
	scissor::Float64
	vcmap::vcMap
	ijmap::ijMap
	band::Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}
	Kernal::KT
end
function Base.show(io::IO, bse::BSEcluster_spinful)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
function BSEcluster_spinful(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction;
	kgrid::RedKgrid = RedKgrid(MonkhorstPack([1, 1, 1])), v, c, scissor::Real = 0, isqgrid::Bool = false)

	vcmap = vcMap(v, c)
	ijmap = ijMap(numorb(TB))

	band = BAND([ReducedCoordinates(0, 0, 0)], TB; vector = true)[1]
	_sum_wave_is_real!(band)

	Kernal(Val(:initialize))

	return BSEcluster_spinful(TB, Float64(scissor), vcmap, ijmap, band, Kernal)
end
function (bse::BSEcluster_spinful)()
	N = length(bse.vcmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, ReducedCoordinates(0, 0, 0))
end
function (bse::BSEcluster_spinful)(q::ReducedCoordinates)
	N = length(bse.vcmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, q)
end
function (bse::BSEcluster_spinful)(H, q::ReducedCoordinates)
	return _BSE_Hamiltonian!(bse, H)
end
function _BSE_Hamiltonian!(bse::BSEcluster_spinful, H)
	vcmap = bse.vcmap
	kernal = bse.Kernal
	band = bse.band
	scissor = bse.scissor

	N = length(vcmap)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 2:N, i in 1:j-1
		n += 1
		tasks[n] = Threads.@spawn begin
			(v′, c′) = vcmap[i]
			(v, c) = vcmap[j]

			ψc′ = band.vectors[:, c′]
			ψc = band.vectors[:, c]
			ψv′ = band.vectors[:, v′]
			ψv = band.vectors[:, v]

			Kᵈ, Kˣ = kernal(ψc′, ψv′, ψv, ψc)

			H[i, j] = Kᵈ + Kˣ
		end
	end
	for i in 1:N
		n += 1
		tasks[n] = Threads.@spawn begin
			(v, c) = vcmap[i]

			ψc = band.vectors[:, c]
			ψv = band.vectors[:, v]

			Kᵈ, Kˣ = kernal(ψc, ψv, copy(ψv), copy(ψc))

			Δϵ = band.values[c] - band.values[v] + scissor
			H[i, i] = Δϵ + Kᵈ + Kˣ
		end
	end
	wait.(tasks)

	return Hermitian(H, :U)
end
