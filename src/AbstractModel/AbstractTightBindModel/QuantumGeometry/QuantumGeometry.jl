
export QuantumGeometry, QuantumMetric

include("./abelian.jl")
include("./nonabelian.jl")

function QuantumGeometry(kgrid::MonkhorstPack, paras...)
	return QuantumGeometry(RedKgrid(kgrid), paras...)
end
function QuantumGeometry(kgrid::AbstractBrillouinZone, paras...)
	return QuantumGeometry(kgrid.kdirect, paras...)
end
function QuantumGeometry(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel, n = 1)
	return QuantumGeometry(kpoints, TB.H, TB.lattice, TB.orb_location, n)
end

function QuantumMetric(sym::Symbol, paras...; kwargs...)
	return QuantumMetric(Val{sym}, paras...; kwargs...)
end
function QuantumMetric(::Union{Val{:QuantumGeometry}, Val{:QG}}, paras...)
	QG = QuantumGeometry(paras...)
	return QuantumMetric(QG)
end
function QuantumMetric(QG::AbstractArray{<:Number, 3})
	return real.(conj.(QG) .+ QG) ./ 2
end
function QuantumMetric(QG::AbstractArray{<:Number, 5})
	#g is the symmetric part of QG
	g = similar(QG)
	for k in axes(QG, 5), ν in 1:3, μ in 1:3
		T = QG[:, :, μ, ν, k]
		g[:, :, μ, ν, k] = (T .+ T') ./ 2
	end
	return g
end


