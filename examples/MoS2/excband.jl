
using JLD2
using KrylovKit
using WannierExcitonModel

period = [1, 1, 0]
unitfolder = "./unitcell"

cell = read(joinpath(unitfolder, "POSCAR"), POSCAR; period)
# Mo:1:10, S1:11:16, S2:17:22
atom_orb_index = [
	collect(1:10),
	collect(11:16),
	collect(17:22),
]
# contruct a set of Wannier pairs.
pair = atom2orbpair(cell, 2, atom_orb_index)

nk = 36
kgrid = MonkhorstPack([nk, nk, 1])

Kᵈ_U = load(joinpath(unitfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(unitfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(unitfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(unitfolder, "pair_bare.jld2"))["interaction"]
# here `pair` means that we only use a subset of interaction to calculate BSE.
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair; pair)
Kˣ = PairInterAction(Kˣ_U, Kˣ_pair; pair)

TB = TightBindModel(; cell,
	hops = joinpath(unitfolder, "wannier90_hr.dat"),
	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
)

# DFTgap:1.6eV; GWgap: 2.8eV
# make sure setting channel=1 for noncollinear case.
bse = BSE(kgrid, TB, Kᵈ, Kˣ; channel=1, v = 13:14, c = 15:16, scissor = 1.2)
kline = Kline(30, bse)

eigsolveconfigure!(; verbosity = 0, howmany = 80, which = EigSorter(abs))
band = BAND(kline, bse)
write("./excband.dat", band.band, BAND_dat; kline)
