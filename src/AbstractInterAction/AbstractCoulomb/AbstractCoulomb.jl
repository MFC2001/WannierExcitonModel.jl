
"""
    abstract type AbstractCoulomb <: AbstractInterAction end

Abstract type for Coulomb interaction potentials.
"""
abstract type AbstractCoulomb <: AbstractInterAction end
"""
    abstract type RealCoulomb <: AbstractCoulomb end

The Coulomb potential defined in real space.
"""
abstract type RealCoulomb <: AbstractCoulomb end
"""
    abstract type ReciprocalCoulomb <: AbstractCoulomb end

The Coulomb potential defined in reciprocal space.
"""
abstract type ReciprocalCoulomb <: AbstractCoulomb end

#1e-19
const qₑ = 1.602176634
#1e-12
const ϵ₀ = 8.854187817
#This term equal to e^2/4πϵ₀ with the unit of r is Å.
#Coulomb potential is CoulombScale/r, the unit of r is Å, potential Energy unit is eV.
const CoulombScale = qₑ * 1e3 / (4 * π * ϵ₀)


include("./RealCoulomb.jl")
include("./ReciprocalCoulomb.jl")
