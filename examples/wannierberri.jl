# work as a post-processing of wannierberri
# generate date file with the format of Wannier90.


using WannierExcitonModel

period = [1, 1, 1]
datafolder = "./"

# read file
xmldata = read(joinpath(datafolder, "pwscf.save/data-file-schema.xml"), QE_xml_new; period)
chkdata = read(joinpath(datafolder, "w90data.chk.npz"), wannierberri_chk)
bands = read(joinpath(datafolder, "w90data.eig.npz"), wannierberri_eig)

# construct k-grid used in MLWF process.
kgrid = RedKgrid(MonkhorstPack(chkdata.mp_grid))
kgrid.kdirect .= xmldata["kpoints"]

cell = xmldata["cell"]
lattice = cell.lattice

# evaluate hopping terms in real space.
hr = HK2HR(build_wannier_hk(bands, chkdata), kgrid, lattice; grid_atol=1e-6)
write(joinpath(datafolder, "wannier90_hr.dat"), hr, wannier90_hr)

# generate wannier90_centres.xyz
orbital = wannier90_centres(map(x -> lattice * x, chkdata.wannier_centres);
    name=fill("X", chkdata.num_wann),
    atom_location=map(x -> lattice * x, cell.location),
    atom_name=cell.name,
)
write(joinpath(datafolder, "wannier90_centres.xyz"), orbital, wannier90_centres)

# calculate band
TB = TightBindModel(; cell, hops=hr, orbital)
# modify kline according to individual needs
kline = [
    "L" => [0.5, 0.5, 0.5],
    "Γ" => [0, 0, 0],
    "X" => [0.5, 0, 0.5],
    "W" => [0.5, 0.25, 0.75],
    "L" => [0.5, 0.5, 0.5],
]
kline = Kline(kline, 50, TB)
band = BAND(kline, TB)
write("./eleband.dat", band, BAND_dat; kline)
