module WannierExcitonModel

using WannierInterActionBase

using WannierInterActionBase.Requires

using WannierInterActionBase.SpecialFunctions

using WannierInterActionBase.StaticArrays
using WannierInterActionBase.StructEquality

Base_exports = names(WannierInterActionBase; all = false)
for sym in Base_exports
    @eval export $sym
end

using LinearAlgebra

using Printf
using Serialization


include("./Topology/Topology.jl")
abstract type AbstractModel end
include("./AbstractTightBindModel/AbstractTightBindModel.jl")
include("./AbstractBSE/AbstractBSE.jl")
include("./HR/HR.jl")
include("./IdealLattice/IdealLattice.jl")
include("./BAND/BAND.jl")
include("./Exciton/Exciton.jl")
# include("./Spglib/Spglib.jl")
include("./Symmetry/Symmetry.jl")
include("./WaveFunction/WaveFunction.jl")
include("./Tools/Tools.jl")



function __init__()
	@require KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77" include("../ext/KrylovKit/KrylovKit.jl")
	@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/Plots/Plots.jl")
end

end
