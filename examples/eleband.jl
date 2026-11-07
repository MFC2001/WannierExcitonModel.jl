# calculate Wannier-interpolating band

using WannierExcitonModel

TB = IdealTB(:hexagonal; a = 2.5, c = 20, t₁ = 1)
# unitfolder = "./unitcell/"
# TB = TightBindModel(;
# 	hops = joinpath(unitfolder, "wannier90_hr.dat"),
# 	cell = joinpath(unitfolder, "POSCAR"),
# 	orbital = joinpath(unitfolder, "wannier90_centres.xyz"),
# 	period = [1, 1, 0],
# )

kline = Kline(30, TB)
# kline = [
# 	"L" => [0.5, 0.5, 0.5],
# 	"G" => [0, 0, 0],
# 	"X" => [0.5, 0, 0.5],
# 	"W" => [0.5, 0.25, 0.75],
# 	"L" => [0.5, 0.5, 0.5],
# ]
# kline = Kline(kline, 30, TB)

band = BAND(kline, TB)

# plot band
using Plots
plotlyjs()
plot(kline, band)
