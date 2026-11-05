
using JLD2
using WannierExcitonModel

period = [1, 1, 0]
unitfolder = "./unitcell"

cell = read(joinpath(unitfolder, "POSCAR"), POSCAR; period)
TB_up = TightBindModel(; cell,
	hops = joinpath(unitfolder, "up", "wannier90_hr.dat"),
	orbital = joinpath(unitfolder, "up", "wannier90_centres.xyz"),
)
TB_dn = TightBindModel(; cell,
	hops = joinpath(unitfolder, "dn", "wannier90_hr.dat"),
	orbital = joinpath(unitfolder, "dn", "wannier90_centres.xyz"),
)

nk = 19
kgrid = MonkhorstPack([nk, nk, 1])

Kᵈ_U = load(joinpath(unitfolder, "U_1_mirror.jld2"))["Uwithlr"]
Kᵈ_pair = load(joinpath(unitfolder, "pair_1.jld2"))["interaction"]
Kᵈ = PairInterAction(Kᵈ_U, Kᵈ_pair)
Kˣ = PairInterAction(2) # This is an empty interaction, because magnon doesn't need Kˣ. nw = nw_up+nw_dn = 2

bse = BSE(kgrid, TB_up, TB_dn, Kᵈ, Kˣ; v_up=[1], c_up=Int[], v_dn=Int[], c_dn=[1])
kline = Kline(30, bse)
# for v_dn to c_up, the total spin of excitation is -1
band = BAND(kline, bse).bandSzn1

# ΔE_gap_g0w0 ≈ 3
# minimum(band) ≈ -3.38
band .-= minimum(band)
write("./bseband.dat", band, BAND_dat; kline)
