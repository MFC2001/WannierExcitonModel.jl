
include("./WilsonLoop.jl")

function BerryCurvature(sym::Symbol, args...; kwargs...)
	return BerryCurvature(Val(sym), args...; kwargs...)
end

function BerryCurvature(::Union{Val{:QuantumGeometry}, Val{:QG}}, args...)
	QG = QuantumGeometry(args...)
	return BerryCurvature(QG)
end
function BerryCurvature(QG::AbstractArray{<:Number, 3})
	return imag.(conj.(QG) .- QG)
end
function BerryCurvature(QG::AbstractArray{<:Number, 5})
	#F is the antisymmetric part of QG.
	F = similar(QG)
	for k in axes(QG, 5), ν in 1:3, μ in 1:3
		T = QG[:, :, μ, ν, k]
		F[:, :, μ, ν, k] = 1im * (T .- T')
	end
	return F
end
