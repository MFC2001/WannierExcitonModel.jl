
using KrylovKit
using WannierExcitonModel

unitfolder = "./unitcell"
TB = TightBindModel(;
	hops = joinpath(unitfolder, "wannier90_hr.dat"),
	cell = joinpath(unitfolder, "POSCAR"),
	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
	period = [1, 1, 1],
)

kgrid = MonkhorstPack([10, 10, 10]; kshift=[1 // 2, 1 // 2, 1 // 2])

Kᵈ_U = load(joinpath(unitfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kˣ_U = load(joinpath(unitfolder, "U_bare_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(unitfolder, "pair_1.jld2"))["interaction"]
Kˣ_pair = load(joinpath(unitfolder, "pair_bare.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair)
Kˣ = PairInterAction(Kˣ_U, Kˣ_pair)
#GW gap ≈ 14.33
bse = BSE(kgrid, TB, Kᵈ, Kˣ; v = [1, 2, 3], c = [4])

kline = [
    "L" => [0.5, 0.5, 0.5],
    "G" => [0, 0, 0],
    "X" => [0.5, 0, 0.5],
    "W" => [0.5, 0.25, 0.75],
    "L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 30, bse)

eigsolveconfigure!(; verbosity=0, howmany=50, which=EigSorter(abs))
band = BAND(kline, bse)

write("./singlet.dat", band.bands, BAND_dat; kline)
write("./triplet.dat", band.bandt, BAND_dat; kline)
