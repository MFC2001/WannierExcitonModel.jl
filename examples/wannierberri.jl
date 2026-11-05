
using WannierExcitonModel

period = [1, 1, 1]
datafolder = "./"

xmldata = read(joinpath(datafolder, "pwscf.save/data-file-schema.xml"), QE_xml_new; period)
chkdata = read(joinpath(datafolder, "w90data.chk.npz"), wannierberri_chk)
bands = read(joinpath(datafolder, "w90data.eig.npz"), wannierberri_eig)

kgrid = RedKgrid(MonkhorstPack(chkdata.mp_grid))
kgrid.kdirect .= xmldata["kpoints"]

hr = HK2HR(build_wannier_hk(bands, chkdata), kgrid, xmldata["cell"].lattice; grid_atol=1e-6)
write(joinpath(datafolder, "wannier90_hr.dat"), hr, wannier90_hr)

cell = convert(xmldata["cell"])
orbital = wannier90_centres(map(x -> cell.lattice * x, chkdata.wannier_centres);
    name=fill("X", chkdata.num_wann),
    atom_location=cell.location,
    atom_name=cell.name,
)
write(joinpath(datafolder, "wannier90_centres.xyz"), orbital, wannier90_centres)


TB = TightBindModel(; cell, hops=hr, orbital)

kline = [
    "L" => [0.5, 0.5, 0.5],
    "Γ" => [0, 0, 0],
    "X" => [0.5, 0, 0.5],
    "W" => [0.5, 0.25, 0.75],
    "L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 50, TB)

band = BAND(kline, TB; vector=true)
write("./eleband.dat", band, BAND_dat; kline)

# using Plots
# plotlyjs()
# plot(kline, band)
