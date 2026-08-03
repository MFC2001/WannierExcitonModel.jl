
using JLD2
using KrylovKit
using WannierExcitonModel

unitfolder = "./unitcell"
TB = TightBindModel(;
	hops = joinpath(unitfolder, "wannier90_hr.dat"),
	cell = joinpath(unitfolder, "POSCAR"),
	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
	period = [1, 1, 1],
)

nk = 9
kgrid = MonkhorstPack([nk, nk, nk])
atol = 1e-4
outfolder = "./rundata"

Kᵈ_U = load(joinpath(unitfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(unitfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(unitfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(unitfolder, "pair_bare.jld2"))["interaction"]

kernel = Kernel(kgrid, Kᵈ_U, Kˣ_U; Kᵈ_pair, Kˣ_pair, atol)
bse = BSE(:SU2, TB, kernel; v = [1, 2, 3], c = [4], scissor = 0)
eigsolveconfigure!(; verbosity = 0, howmany = 50, which = EigSorter(abs))

nk = 41
kgrid = MonkhorstPack([nk, nk, 1])
BC = BerryCurvature(Val(:WilsonLoop), kgrid, TB.lattice)

(uband_t, uband_s) = BAND(BC.vertex, bse; vector = true, wfctype = :Periodic)

uband_t_0 = map(x -> x(0, [1, 2, 3]), uband_t)
# uband_t_1 = map(x -> x(1, [1, 2, 3]), uband_t)

u_getter = WilsonLoopWaveGetter(map(U -> U.vectors[:, [1]], uband_t_0))
BC(:CN, u_getter)

using Plots
plotlyjs()

plot(BC)
# plot(BC; colorbar = true)
