
nk = 41
kgrid = MonkhorstPack([nk, nk, 1])
BC = BerryCurvature(Val(:WilsonLoop), kgrid, bse.lattice)

band = BAND(BC.vertex, bse; vector = true, wfctype = :Periodic)

# 0 is the weight of exciton center coordinate, can be in [0, 1]
# [1, 2, 3] is the band index
# uband_0 is a `Vector` of `Eigen`.
uband_0 = map(x -> x(0, [1, 2, 3]), band.uband)
# uband_1 = map(x -> x(1, [1, 2, 3]), band.uband)

u_getter = WilsonLoopWaveGetter(map(U -> U.vectors[:, [1]], uband_0))
BC(:CN, u_getter)

using Plots
plotlyjs()
plot(BC)
