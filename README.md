This repository contains the code and all raw- and metadata used in the the manuscript:
<b> Buchwald SZ, Frank AB, Birgel D, Senger K, Mosociova T, Pei Y, Beaty B, Tarhan LG, Galasso F, Gómez Correa MA, Grasby SE, Struck U,  Steinkrauss R, Gliwa J, Lahajnar N, Peckmann J, and Foster WJ
<br> "Reconstructing environmental and microbial ecosystem changes across the Permian–Triassic mass extinction at Lusitaniadalen, Svalbard "

</B> To properly run the script, the folder structure must mirror the structure in this repository:
The R scripts `Lusitaniadalen.R` and `Lusitaniadalen_Functions.R` need to be in the working directory.
Create a sub-folder titled `Raw_data` in the working directory, which should contain all raw- and metadata, including
`d13Corg.txt`, `Geochem.txt`, `HAWK.txt`, `Lusitaniadalen_Metadata.xlsx`, and the folders `area_molecular_sieve` and `area_n_alk`, in which the integrated
peak areas of the compounds of interest are saved in an individual `.txt`-files per sample.

When initiating the project, a folder `Output` is created in the working directory, that will contain all output produced.

Additionally, the output of the calculation of element enrichment factors is provided in the Excel file `geochem_out.xlsx`, quantified n-alkanes (in ug/g TOC) in the Excel file `FID_n_alk_quantified_ug_TOC_Chol.xlsx`, and quantified isoprenoids (in ug/g TOC) in the Excel file `FID_mol_quantified_ug_TOC_Chol.xlsx`. 
