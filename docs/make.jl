using Documenter, WannierExcitonModel

pages = [
	"Home" => "index.md",
	"Getting Started" => "start.md",
	"Model" => [
		"Tight Binding Model" => "tbmodel.md",
		"Excitonic BSE Model" => "bsemodel.md",
	],
	"Structure" => "structure.md",
	"Reciprocal Hopping" => "rh.md",
	"Brillouin Zone" => "brillouin_zone.md",
	"Interaction" => [
		"Coulomb" => "coulomb.md",
		"Mirror Correction" => "mirror_correction.md",
		"UwithLR" => "uwithlr.md",
		"Kernal" => "kernal.md",
	],
	"Topology" => "topology.md",
	"IO" => "io.md",
]

format = Documenter.HTML(;
	prettyurls = true,
	collapselevel = 1,
)

makedocs(;
	# modules = [WannierExcitonModel],
	sitename = "WannierExcitonModel.jl",
	authors = "Fanchen Meng",
	format,
	pages,
)

deploydocs(
    repo = "github.com/MFC2001/WannierExcitonModel.jl.git",
)
