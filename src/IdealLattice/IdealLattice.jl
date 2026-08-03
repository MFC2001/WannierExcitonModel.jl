export IdealTB

"""
    IdealTB(sym::Symbol, args...; kwargs...)

See details in [`IdealCell`](@ref) and [`IdealHR`](@ref).
"""
IdealTB(sym::Symbol, args...; kwargs...) = IdealTB(Val(sym), args...; kwargs...)

include("./GammaGraphyne.jl")
include("./hexagonal.jl")
include("./kagome_hexagonal.jl")
include("./kagome.jl")
include("./star.jl")
