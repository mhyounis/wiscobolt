These elastic scattering cross sections were produced by the ELSEPA program, version 2020,
created by Francesc Salvat, Aleksander Jablonski, and Cedric J. Powell.

ELSEPA can be downloaded at : https://doi.org/10.17632/w4hm5vymym.1

A copy of ELSEPA's Apache license is included ("LICENSE.txt").

None of the data files produced by ELSEPA have been modified in any way. ELSEPA input files corresponding to
the the elements are all provided in the element folders, so that the user may understand what settings generated
these cross sections.

ELSEPA calculates elastic scattering cross sections for incident relativistic electrons or positrons
upon atoms, molecules, and ions. A brief description of the method is given here, however, the user is referred
to the following three sources:
	[1] Salvat et al. elsepa—Dirac partial-wave calculation of elastic scattering of electrons and positrons by atoms, positive ions and molecules. (2005)
	
	[2] Salvat. PENELOPE-2018: A code system for Monte Carlo simulation of electron and poton transport. (2019)

	wiscobolt physics (included in wiscobolt/Documents folder, describes how wiscobolt discretizes the cross sections provided)

Briefly, ELSEPA relies on a program from the same authors, RADIAL, that solves the radial Dirac equation for atomic number Z.
Then, partial wave analysis is performed in order to construct differential scattering cross sections for a chosen projectile.
In order to treat molecules, the independent atom approximation is used, wherein molecular bond structure is neglected
and thus molecules are considered to be a grouping of independent atoms which are rigidly separated from one another.
This results in a molecular cross section which is a sum of independent atom approximations but also a sum over every
pair of atoms of the product of their complex scattering amplitudes and another factor relating to momentum transfer
and the distance between atoms.

The cross sections are given between electron energies of 10 eV and 100 MeV (100,000,000 eV). There are 73 energies.
Increment sizes are:
Energy range				Energy increments
10 eV - 100 eV: 			5 eV
100 eV - 1,000 eV: 			100 eV
1,000 eV - 10,000 eV: 			1,000 eV
10,000 eV - 100,000 eV: 		10,000 eV
100,000 eV - 1,000,000 eV: 		100,000 eV
1,000,000 eV - 10,000,000 eV: 		1,000,000 eV
10,0000,00 eV - 100,000,000 eV: 	10,000,000 eV
The specific energies are listed in "Evals.txt"

The data files are named as "dcs_WpXYZe0A.dat", which designates the cross section for
electron energy W.XYZ*10^A eV.