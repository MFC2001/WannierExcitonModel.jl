
"""
    abstract type AbstractInterAction end

Abstract type for interaction potentials. 
"""
abstract type AbstractInterAction end

include("./AbstractCoulomb/AbstractCoulomb.jl")
include("./LRCorrection/LRCorrection.jl")

include("./MirrorCorrection/MirrorCorrection.jl")
