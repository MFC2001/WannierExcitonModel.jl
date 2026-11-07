# calculate Berry curvature and Chern Number

using WannierExcitonModel

TB = IdealTB(:hexagonal; a = 2.5, c = 20, t₁ = 1, tₕ = -0.1)
# unitfolder = "./unitcell/"
# TB = TightBindModel(;
# 	hops = joinpath(unitfolder, "wannier90_hr.dat"),
# 	cell = joinpath(unitfolder, "POSCAR"),
# 	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
# 	period = [1, 1, 0],
# )

# construct a `BerryCurvature`, which contains all the small elements and the closed paths used for calculating the Wilsonloop.
nk = 41
kgrid = MonkhorstPack([nk, nk, 1])
BC = BerryCurvature(Val(:WilsonLoop), kgrid, TB.lattice)

# wfctype = :Periodic means using the atomic guage.
uband = BAND(BC.vertex, TB; vector = true, wfctype = :Periodic)

# `WilsonLoopWaveGetter` can provide wavefunctions to calculate Wilsonloop
# [1] means calculating the band of index 1
# it supports multi-band case
u_getter = WilsonLoopWaveGetter(map(U -> U.vectors[:, [1]], uband))
# Chern Number
chern = BC(:CN, u_getter)

# plot the Berry curvature in reciprocal space
using Plots
plotlyjs()
plot(BC)
