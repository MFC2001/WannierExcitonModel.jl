export AbstractBSE
"""
    AbstractBSE <: AbstractModel

Abstract type used to construct excitonic BSE Model.
"""
abstract type AbstractBSE <: AbstractModel end

include("./BSEModel/BSEModel.jl")
include("./_BSE_berryconnection_vck.jl")
include("./shell.jl")
include("./vck2ijR/vck2ijR.jl")
include("./BSE_uEigen/BSE_uEigen.jl")

include("./others/others.jl")
