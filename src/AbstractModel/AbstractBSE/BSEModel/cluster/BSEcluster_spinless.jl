struct BSEcluster_spinless{
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
function Base.show(io::IO, bse::BSEcluster_spinless)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
function BSEcluster_spinless(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction;
	kgrid::RedKgrid = RedKgrid(MonkhorstPack([1, 1, 1])), v, c, scissor::Real = 0, isqgrid::Bool = false)

	vcmap = vcMap(v, c)
	ijmap = ijMap(numorb(TB))

	band = BAND([ReducedCoordinates(0, 0, 0)], TB; vector = true)[1]
	_sum_wave_is_real!(band)

	Kernal(Val(:initialize))

	return BSEcluster_spinless(TB, Float64(scissor), vcmap, ijmap, band, Kernal)
end
function (bse::BSEcluster_spinless)()
	N = length(bse.vcmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)
	return bse(Htriplet, Hsinglet, ReducedCoordinates(0, 0, 0))
end
function (bse::BSEcluster_spinless)(q::ReducedCoordinates)
	N = length(bse.vcmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)
	return bse(Htriplet, Hsinglet, q)
end
function (bse::BSEcluster_spinless)(Htriplet, Hsinglet, q::ReducedCoordinates)
	return _BSE_Hamiltonian!(bse, Htriplet, Hsinglet)
end
function _BSE_Hamiltonian!(bse::BSEcluster_spinless, Htriplet, Hsinglet)
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

			Htriplet[i, j] = Kᵈ
			Hsinglet[i, j] = Kᵈ + 2 * Kˣ
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
			Htriplet[i, i] = Δϵ + Kᵈ
			Hsinglet[i, i] = Δϵ + Kᵈ + 2 * Kˣ
		end
	end
	wait.(tasks)

	return Hermitian(Htriplet, :U), Hermitian(Hsinglet, :U)
end
