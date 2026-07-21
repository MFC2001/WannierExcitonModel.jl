export NP2SP
function NP2SP(hr::HR, orbital::wannier90_centres;  buildindex = true)

	spinhr = spinHR(hr; buildindex)
	spinorbital = spin_centres(orbital)

	return spinhr, spinorbital
end

