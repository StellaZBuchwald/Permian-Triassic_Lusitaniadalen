# This script performs all calculations from the manuscript:
# Buchwald SZ, Frank AB, Birgel D, Senger K, Mosociova T, Pei Y, Tarhan L, Galasso F, Gómez Correa MA, 
# Grasby SE, Struck U, Steinkrauss R, Gliwa J,  Lahajnar N, Peckmann J, and Foster WJ

# "Reconstructing environmental and microbial ecosystem changes across the Permian–Triassic mass
# extinction at Lusitaniadalen, Svalbard"

# This script and the script "Lusitaniadalen_Functions.R" need to be in the working directory to be properly
# loaded. Create a sub-folder titled "Raw_data" in the working directory, which should contain all raw 
# and metadata, including "Geochem.txt", "d13Corg.txt", "HAWK.txt", "Homohopanes.txt",
# Lusitaniadalen_Metadata.xlsx", and the folders "area_mol_sieve" and "area_n_alk", in which the
# integrated peak areas of the compounds of interest are saved in an individual .txt-files per sample.

# When initiating the project, a folder "Output" is created in the working directory,
# that will contain all data produced.

setwd("") # set working directory
source('Lusitaniadalen_Functions.R')
init_project()

#_______________________________________________________________________________
#                          (1) GEOCHEMICAL DATA
#_______________________________________________________________________________
# calculate redox-sensitive metal enrichments as well as standard-normalized
# rare Earth element data including Ce anomalies.

# import data
geochem <- read.table("Raw_data/Geochem.txt", sep = "\t", header = TRUE)

# Al-normalized concentrations for V, Mo, URe 
geochem$VAl <- geochem$V / geochem$Al
geochem$MoAl <- geochem$Mo / geochem$Al 
geochem$UAl <- geochem$U / geochem$Al
geochem$ReAl <- geochem$Re / geochem$Al

# PAAS-normalized Enrichment Factors for V, Mo, U and Re
Al_PAAS <- 10
geochem$VEF <- geochem$VAl / (140/Al_PAAS)
geochem$MoEF <- geochem$MoAl / (3/Al_PAAS)
geochem$UEF <- geochem$UAl / (3.1/Al_PAAS)
geochem$ReEF <- geochem$ReAl / (0.0004/Al_PAAS)

# adjusted PAAS-normalized Enrichment Factors for V, Mo, U and Re
Al_PAAS <- 10
geochem$eV <- geochem$V - (geochem$Al * (140/Al_PAAS))
geochem$eMo <- geochem$Mo - (geochem$Al * (3/Al_PAAS))
geochem$eU <- geochem$U - (geochem$Al * (3.1/Al_PAAS))
geochem$eRe <- geochem$Re - (geochem$Al * (0.0004/Al_PAAS))

geochem$cVEF <- (geochem$eV + 140)/140
geochem$cMoEF <- (geochem$eMo + 3)/3
geochem$cUEF <- (geochem$eU + 3.1)/3.1
geochem$cReEF <- (geochem$Re + 0.0004)/0.0004

# PAAS-normalized REE+Y, Ce & Y anomalies and LREE 
geochem$LaPAAS <- geochem$La / 38
geochem$CePAAS <- geochem$Ce / 80
geochem$PrPAAS <- geochem$Pr / 8.9
geochem$NdPAAS <- geochem$Nd / 32
geochem$SmPAAS <- geochem$Sm / 5.6
geochem$EuPAAS <- geochem$Eu / 1.1
geochem$GdPAAS <- geochem$Gd / 4.7
geochem$TbPAAS <- geochem$Tb / 0.77
geochem$DyPAAS <- geochem$Dy / 4.4
geochem$YPAAS <- geochem$Y / 27
geochem$HoPAAS <- geochem$Ho / 1
geochem$ErPAAS <- geochem$Er / 2.9
geochem$TmPAAS <- geochem$Tm / 0.4
geochem$YbPAAS <- geochem$Yb / 2.8
geochem$LuPAAS <- geochem$Lu / 0.43
geochem$CeCe <- geochem$CePAAS / (geochem$PrPAAS * (geochem$PrPAAS / geochem$NdPAAS))
geochem$YHo <- geochem$YPAAS / geochem$HoPAAS
geochem$LREE <- geochem$SmPAAS / geochem$YbPAAS

write.xlsx(geochem, file = "Output/geochem_out.xlsx")

#_______________________________________________________________________________
#                     (2) ORGANIC GEOCHEMISTRY
#_______________________________________________________________________________

############################### DATA IMPORT ####################################

# quantify n-alkanes (ug compound/g TOC)
a_nalk <- read_raw_data("Raw_data/area_n_alk")
q_nalk_ug_TOC <- quantify_compounds(a_nalk, "meta_nalk")

# quantify post molecular sieve comppunds
a_mol <- read_raw_data("Raw_data/area_molecular_sieve")
q_mol_ug_TOC <- quantify_compounds(a_mol, "meta_mol")

# homohopanes
# peak areas of m/z = 191 from molecular sieve treated samples
homohopanes <- read.table("Raw_data/Homohopanes.txt", sep = "\t", header = TRUE)
homohopanes$HHI <- homohopanes$C31S/(homohopanes$C31S + homohopanes$C31R)
mean(homohopanes$HHI)

# metadata
meta_nalk <- read.xlsx("Raw_data/Lusitaniadalen_Metadata.xlsx", sheet = "meta_nalk")
meta_mol <- read.xlsx("Raw_data/Lusitaniadalen_Metadata.xlsx", sheet = "meta_mol")

# TOC data
TOC <- merge(as.data.frame(cbind(log_height = meta_nalk$log_height, TOC_acid = meta_nalk$TOC)),
             as.data.frame(cbind(log_height = read.table("Raw_data/HAWK.txt", sep = "\t", head = TRUE)$log_height,
                                 TOC_HAWK = read.table("Raw_data/HAWK.txt", sep = "\t", head = TRUE)$TOC)),
             by = "log_height")
# compare TOC
TOC$deltaTOC <- TOC$TOC_HAWK - TOC$TOC_acid
mean(TOC$deltaTOC)

# d13C ORG
d13C <- read.table("Raw_data/d13Corg.txt", sep = "\t", header = TRUE)

###################### DATA PREPARATION FOR PLOTTING ###########################

# N-ALKANES
nalk <- merge(q_nalk_ug_TOC, as.data.frame(select(meta_nalk, sample, log_height)), by = "sample") # merge log_height
write.xlsx(nalk, file = "Output/FID_n_alkanes_quantified_ug_TOC_Chol.xlsx")

# Terrigenous-aquatic ratio (TAR)
TAR <- cbind(TAR = (nalk$`n-C27` + nalk$`n-C29` + nalk$`n-C31`)/(nalk$`n-C15` + nalk$`n-C17` + nalk$`n-C19`),
             log_height = nalk$log_height)

# Carbon preference index (CPI)
nalk[is.na(nalk)] <- 0 # replace NAs by 0
CPI <- as.data.frame(cbind(CPI = ((((nalk$`n-C25` + nalk$`n-C27` + nalk$`n-C29` + nalk$`n-C31` + nalk$`n-C33`)/
            (nalk$`n-C24` + nalk$`n-C26` + nalk$`n-C28` + nalk$`n-C30` + nalk$`n-C32`)) +
           ((nalk$`n-C25` + nalk$`n-C27` + nalk$`n-C29` + nalk$`n-C31` + nalk$`n-C33`)/
              (nalk$`n-C26` + nalk$`n-C28` + nalk$`n-C30` + nalk$`n-C32` + nalk$`n-C34`)))/2),
           log_height = nalk$log_height))

# POST MOLECULAR SIEVE COMPOUNDS
mol <- merge(q_mol_ug_TOC, as.data.frame(select(meta_mol, sample, log_height)), by = "sample") # merge log_height
write.xlsx(mol, file = "Output/FID_molsieb_quantified_ug_TOC_Chol.xlsx")

# pristane and phytane
Pr_Ph <- melt(as.data.frame(select(mol, log_height, Pristane, Phytane)), id.vars = 'log_height', variable.name = 'compound')
mol$Pr_Ph_r <- mol$Pristane/mol$Phytane #pristane/phytane

# isoprenoids
mol[is.na(mol)] <- 0 # replace NAs by 0
mol$sumC21_24_isos <- mol$`C21-isoprenoid`+ mol$`C22-isoprenoid` + mol$`C23-isoprenoid` + mol$`C24-isoprenoid`

PMI <- melt(as.data.frame(cbind(log_height = mol$log_height,
                                PMI = mol$PMI_regular,
                                sum_pseudos = mol$sumC21_24_isos)), id.vars = 'log_height', variable.name = 'compound')

short_isos <- melt(as.data.frame(cbind(log_height = mol$log_height,
                                    "C14" = mol$`C14-isoprenoid`,
                                    "C15" = mol$farnesane,
                                    "C16" = mol$`C16-isoprenoid`,
                                    "C17" = mol$`C17-isoprenoid`,
                                    "C18" = mol$`nor-pristane`)),  id.vars = 'log_height', variable.name = 'compound')

################################ PLOTTING ######################################

plot_d13C <- ggplot(d13C, aes(log_height, d13Corg)) +
  geom_point(aes(colour = "d13Corg"), pch = 1, size = 4, stroke = 1, show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "height (m)",
                         y_name = expression(paste(delta^"13", "C"[org], " (\u2030)"))) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90))

############ Carbon Preference Index (CPI)
plot_CPI <- ggplot(CPI, aes(log_height, CPI)) +
  geom_point(aes(colour = "CPI"), pch = 1, size = 4, stroke = 1,  show.legend = FALSE) +
  geom_hline(yintercept = 1.0, colour = "grey13", linetype = "dashed") + # vertical line at CPI = 1
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "CPI") +
  ylim(0.5, 1.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text( angle = 90))

############ Terrigenous-aquatic ratio (TAR)
plot_TAR <- ggplot(TAR, aes(log_height, TAR)) +
  geom_point(aes(colour = "TAR"), pch = 1, size = 4, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = "TAR"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "TAR") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text( angle = 90))

############ Pristane + phytane
plot_PrPh <- ggplot(Pr_Ph, aes(log_height, value)) +
  geom_point(aes(colour = compound, shape = compound), size = 4, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(name="",values = c("olivedrab4", "olivedrab2"))+
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "") +
  scale_shape_manual(values = c(4, 1)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "bottom")

############ Pristane/Phytane
plot_PrPh_ratio <- ggplot(mol, aes(log_height, Pr_Ph_r)) +
  geom_point(pch = 1, size = 4, stroke = 1, show.legend = FALSE) +
  geom_line(linetype = "dotted", show.legend = FALSE) +
  geom_hline(yintercept = 1.0, colour = "grey13", linetype = "dashed") + # vertical line at Pr/Ph = 1
  scale_color_manual(values = c("grey13", "grey70")) +
  plot_common_parameters(y_name = "Pr/Ph",
                         x_name = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text( angle = 90))

# Regular PMI + sum pseudos
plot_PMI <- ggplot(PMI, aes(log_height, value)) +
  geom_point(aes(colour = compound, shape = compound), size = 4, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(name="",values = c("grey13", "grey70"))+
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "")+
  scale_shape_manual(values = c(1, 16)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "bottom")

# shorter isoprenoids
col <- paletteer_d("fishualize::Callanthias_australis") # define color palette
plot_isoprenoids <- ggplot(short_isos, aes(log_height, value)) +
  geom_point(aes(colour = compound), pch = 20, size = 4, stroke = 1, show.legend = FALSE) +
  scale_color_manual(name="",values = col) +
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text(angle = 90),
        legend.position = "right")

############ Total Plots in one figure
total_biomarker <- ggarrange(plot_d13C, plot_CPI, plot_TAR, plot_PrPh_ratio,  plot_PrPh, plot_PMI, plot_isoprenoids,
                  nrow = 1,  align = "hv")
total_biomarker
ggsave(total_biomarker, file = "Output/Total_biomarker.pdf", width = 70, height = 55, units = "cm")
