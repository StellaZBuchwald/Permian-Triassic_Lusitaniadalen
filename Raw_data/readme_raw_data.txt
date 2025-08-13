The Raw data folder contains several files and sub-folder:

(A)
area_nalk
Integrated peak areas of the n-alkanes are saved in an individual .txt-files per sample.

(B)
area_molecular_sieve
Integrated peak areas of the compounds of interest after molecular sieve treatment are
saved in an individual .txt-files per sample.

(C)
Homohopane.txt
Integrated peak areas of R- and S-isomers of the C31-homohopane from m/z = 191 in the
molecular sieve treated samples.

(D)
Pres_Abs_Biphytane_Lycopane_Arylisoprenoids.txt
Contains presence (1)/absence (0) data of the compounds biphytane/lycopane and short-chain
aryl-isoprenoids in the molecular sieve treated samples.

(E)
Lusitaniadalen_Metadata.xlsx
contains two sheets:
(1) "meta-nalk" which contains the metadata for the n-alkanes, and
(2) "meta_mol" which contains the metadata for the molar sieve-treated
samples.
Each sheet contains a table with metadata.

The columns of the sheet "meta_nalk" are:
	sample: sample ID
	log_height: height in the log [m]
	weight: sample weight used for extraction [g]
	TOC: Total organic carbon [wt%]
	a_Chol: area after peak integration with Chromeleon of Cholestane
	c_Chol: concentration of Cholestane used as external standard [ng/uL]
	v_Chol: volume of Cholestane used as external standard [uL]
	v_inject: volumne injected into the FID [uL] 
	v_injec_from: volume from which a certain volume was injected into the FID [uL]
	TLE: fraction of total lipid extract used for maltene-asphaltene separation [%]
	malt: fraction of maltene used for column chromatography [%]
	m_malt: mass of 100% of the maltene fraction [mg]
	F1: fraction of F1 used for addition of external standard [%]
The table "meta_mol" contains the same columns with two additional columns:
	v_istd: volume of the internal standard 10-methyl-nonadecane used [uL]
	c_istd: concentration of the internal standard 10-methyl-nonadecane used [uL][ng/ul]

(F)
d13Corg.txt
Contains stable organic carbon isotope data.
units:	log_height = m
	d13Corg =per mille

(G)
HAWK.txt
Contains all data from HAWK pyrolysis.
units:
	log_height = m
	Tamx = °C
	S0 = mg HC/g rock
	S1 = mg HC/g rock
	S2 = mg HC/g rock
	S3 = mg CO2/g rock
	PC = wt%
	HI = mg HC/g TOC
	OI = mg CO2/g TOC
	TOC = wt%

(H)
geochem.txt
Results from major, minor and trance element analysis as input for calculations of 
redox-sensitive metal enrichments and standard-normalized rare Earth element
data including Ce anomalies. The unit for Al, Ca, Fe, K, Mg, Na, Ti, TOC, P and S is wt%,
the remaining data is given in ppm.
