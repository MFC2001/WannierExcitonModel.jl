
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

kline = [
	"L" => [0.5, 0.5, 0.5],
	"G" => [0, 0, 0],
	"X" => [0.5, 0, 0.5],
	"W" => [0.5, 0.25, 0.75],
	"L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 20, bse)

eigsolveconfigure!(; verbosity = 0, howmany = 50, which = EigSorter(abs))
(band_t, band_s) = BAND(kline, bse; vector = true)

write(joinpath(outfolder, "band_t.dat"), band_t, BAND_dat; kline)
write(joinpath(outfolder, "band_s.dat"), band_s, BAND_dat; kline)
