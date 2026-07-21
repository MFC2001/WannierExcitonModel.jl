struct BSEcluster_general{
	TBT <: AbstractTightBindModel,
	KT <: KernelInterAction,
} <: AbstractBSE
	TB::TBT
	scissor::Float64
	vcmap::vcMap
	ijmap::ijMap
	band::Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}
	kernel::KT
end
function Base.show(io::IO, bse::BSEcluster_general)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
function BSEcluster_general(TB::AbstractTightBindModel, kernel::KernelInterAction;
	kgrid::RedKgrid = RedKgrid(MonkhorstPack([1, 1, 1])), v, c, scissor::Real = 0, isqgrid::Bool = false)

	vcmap = vcMap(v, c)
	ijmap = ijMap(numorb(TB))

	band = BAND([ReducedCoordinates(0, 0, 0)], TB; vector = true)[1]
	phase_normalize!(band, :real_sum)

	kernel(Val(:initialize))

	return BSEcluster_general(TB, Float64(scissor), vcmap, ijmap, band, kernel)
end
function (bse::BSEcluster_general)()
	N = length(bse.vcmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, ReducedCoordinates(0, 0, 0))
end
function (bse::BSEcluster_general)(q::ReducedCoordinates)
	N = length(bse.vcmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, q)
end
function (bse::BSEcluster_general)(H, q::ReducedCoordinates)
	return _BSE_Hamiltonian!(bse, H)
end
function _BSE_Hamiltonian!(bse::BSEcluster_general, H)
	vcmap = bse.vcmap
	kernel = bse.kernel
	band = bse.band
	scissor = bse.scissor

	N = length(vcmap)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 2:N, i in 1:j-1
		n += 1
		tasks[n] = Threads.@spawn begin
			@inbounds begin
				(v′, c′) = vcmap[i]
				(v, c) = vcmap[j]

				ψc′ = band.vectors[:, c′]
				ψc = band.vectors[:, c]
				ψv′ = band.vectors[:, v′]
				ψv = band.vectors[:, v]

				Kᵈ, Kˣ = kernel(ψc′, ψv′, ψv, ψc)
			end

			if isfinite(Kᵈ) && isfinite(Kˣ)
				@inbounds begin
					H[i, j] = Kᵈ + Kˣ
				end
			else
				throw(DomainError(
					(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ, Kˣ = Kˣ),
					"BSE matrix element (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
				))
			end
		end
	end
	for i in 1:N
		n += 1
		tasks[n] = Threads.@spawn begin
			@inbounds begin
				(v, c) = vcmap[i]

				ψc = band.vectors[:, c]
				ψv = band.vectors[:, v]

				Kᵈ, Kˣ = kernel(ψc, ψv, copy(ψv), copy(ψc))
			end

			if isfinite(Kᵈ) && isfinite(Kˣ)
				@inbounds begin
					Δϵ = band.values[c] - band.values[v] + scissor
					H[i, i] = Δϵ + Kᵈ + Kˣ
				end
			else
				throw(DomainError(
					(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ, Kˣ = Kˣ),
					"BSE matrix element (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
				))
			end
		end
	end

	try
		wait.(tasks)
	catch e
		@error "One or more BSE matrix element tasks failed at q=$q" exception = e
		fill!(H, zero(ComplexF64))
		H[diagind(H)] .= ComplexF64(-100.0)
	end

	return Hermitian(H, :U)
end
