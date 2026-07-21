
abstract type _spinaware_Kᵈ end
abstract type _spinaware_Kˣ end

# 如果包含自旋，那么无法定义合适的J, 物理上不合理，故不考虑该情况。
# 一般，若有自旋还需要J-type相互作用，那么wannier基则需要满足时间反演对称性，此时应使用spinblind情况，仅提供无自旋的相互作用。

include("./_U_Kd.jl")
include("./_U_Kx.jl")
include("./_U_pair_Kd.jl")
include("./_U_pair_Kx.jl")
