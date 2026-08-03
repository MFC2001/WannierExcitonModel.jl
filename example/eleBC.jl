
using WannierExcitonModel

TB = IdealTB(:hexagonal; a = 2.5, c = 20, t₁ = 1, tₕ = -0.1)
# unitfolder = "./unitcell/"
# TB = TightBindModel(;
# 	hops = joinpath(unitfolder, "wannier90_hr.dat"),
# 	cell = joinpath(unitfolder, "POSCAR"),
# 	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
# 	period = [1, 1, 0],
# )

nk = 41
kgrid = MonkhorstPack([nk, nk, 1])
BC = BerryCurvature(Val(:WilsonLoop), kgrid, TB.lattice)

uband = BAND(BC.vertex, TB; vector = true, wfctype = :Periodic)

u_getter = WilsonLoopWaveGetter(map(U -> U.vectors[:, [1]], uband))
BC(:CN, u_getter)

using Plots
plotlyjs()

plot(BC)
