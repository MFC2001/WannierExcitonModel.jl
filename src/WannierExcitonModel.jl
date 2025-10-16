module WannierExcitonModel

using Requires

using SpecialFunctions

using LinearAlgebra

using StaticArrays
using StructEquality

using DelimitedFiles
using Printf
using Dates: Dates

using Serialization

include("./Core/Core.jl")
include("./AbstractInterAction/AbstractInterAction.jl")
include("./Topology/Topology.jl")
include("./AbstractModel/AbstractModel.jl")
include("./HR/HR.jl")
include("./IO/IO.jl")
include("./BrillouinZone/BrillouinZone.jl")
include("./BAND/BAND.jl")
include("./Exciton/Exciton.jl")
include("./IdealLattice/IdealLattice.jl")
include("./shell/shell.jl")
# include("./Spglib/Spglib.jl")
include("./Symmetry/Symmetry.jl")
include("./wannier/wannier.jl")
include("./WaveFunction/WaveFunction.jl")
include("./Tools/Tools.jl")



function __init__()
	@require KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77" include("../ext/KrylovKit/KrylovKit.jl")
end

end
