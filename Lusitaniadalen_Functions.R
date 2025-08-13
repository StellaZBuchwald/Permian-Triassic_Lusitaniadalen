library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(paletteer)
library(ggpubr)
library(readtext)
library(reshape2)

# initiate project
init_project <- function() {
    rm(list = ls())
    dir.create("Output/", recursive = TRUE, showWarnings = FALSE)
}

# import raw data (peaka areas)
read_raw_data <- function(folder_path) {
  raw_files <- list.files(path = folder_path)
  data <- data.frame(matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "compound")))
  for (filename in raw_files) {
    fullfilename <- paste(folder_path, filename, sep = "/")
    sample <- read.table(fullfilename, sep = "\t", header = TRUE)
    data <- merge(data, sample, by = "compound", all = TRUE)
  }
  
  rownames(data) <- data$compound
  data <- data[, !names(data) %in% "compound"]
  data <- as.data.frame(t(data))
  
  return(data)
}

# Quantify Compounds
#
# This function takes in a data frame of compound areas and a sheet name from a metadata file.
# It performs quantification calculations on the compounds based on the provided metadata.
# The resulting data frame contains the quantities of compounds normalized to grams of total organic carbon (TOC).
#
# @param a_data A data frame of compound areas.
# @param sheet A character string specifying the sheet name in the metadata file.
#
# @return A data frame with the quantities of compounds normalized to grams of TOC.
#
# @examples
# # Load the metadata file
# meta_nalk_raw <- read.xlsx("Lusitaniadalen_Metadata.xlsx", sheet = "Sheet1")
#
# # Quantify compounds using the function
# quantified_data <- quantify_compounds(a_data, "Sheet1")
#
quantify_compounds <- function(a_data, sheet) {
    # sample: sample ID
    # log_height: height in the log [m]
    # weight: sample weight used for extraction [g]
    # TOC: Total organic carbon [wt%]
    # a_Chol: area after peak integration with Chromeleon of Cholestane
    # c_Chol: concentration of Cholestane used as external standard [pg/uL]
    # v_Chol: volume of Cholestane used as external standard [uL]
    # v_inject: volumne injected into the FID [uL]
    # v_injec_from: volume from which a certain volume was injected into the FID [uL]
    # TLE: fraction of total lipid extract used for maltene-asphaltene separation [%]
    # malt: fraction of maltene used for column chromatography [%]
    # F1: fraction of F1 used for addition of external standard [%]

    meta_nalk_raw <- read.xlsx("Raw_data/Lusitaniadalen_Metadata.xlsx", sheet = sheet) # import metadata

    meta <- meta_nalk_raw[order(match(meta_nalk_raw$sample, rownames(a_data))), ] # order the metadata to the same order as the a_nalk

    # get variables for quantification in vectors
    g_TOC <- (meta$TOC * meta$weight) / 100 # calculate g TOC per g used sediment
    ng_Chol <- meta$c_Chol * meta$v_Chol ## calculate mass of cholestane
    a_Chol <- a_data$Cholestane # area Cholestane (external standard)

    # new data frame with quantities of compounds
    mat <- matrix(NA, nrow = nrow(a_data), ncol = ncol(a_data)) # empty matrix with dimensions of FID_area
    q_data <- data.frame(mat) # matrix to data frame


    # quantification:
    # 1. mass compound = (area_compound*ug_Chol)/area_Chol
    # unit: ng (but calculated only for the fraction of the sample that was injected to the FID)
    for (col in c(1:ncol(a_data))) {
        q_data[, col] <- (a_data[, col] * ng_Chol) / a_Chol
    }

    names(q_data) <- names(a_data)
    rownames(q_data) <- rownames(a_data)

    # 2. mass compound in 100% sample = mass comp * volume during injection [1] * (100/F1 [%]) * (100/malt [%]) * (100/TLE [%])
    # [1] volume during injection = (1/injection volume) * volume from which it was injected
    # unit: ng
    q_data_ng <- q_data

    for (col in c(1:ncol(q_data))) {
        q_data_ng[, col] <- q_data[, col] * (100 / meta$F1) * (100 / meta$malt) * (100 / meta$TLE)
    }

    q_data_ug <- q_data_ng / 1000 # 3. change unit from ng to ug


    # 4. normalize to g TOC
    q_data_ug_TOC <- q_data_ug / g_TOC


    # formatting for export
    q_data_ug_TOC$sample <- rownames(q_data_ug_TOC)
    ncol(q_data_ug_TOC)
    q_data_ug_TOC <- q_data_ug_TOC[, c(ncol(q_data_ug_TOC), 1:ncol(q_data_ug_TOC) - 1)]

    return(q_data_ug_TOC)
}


# Function for customizing plots:
plot_common_parameters <- function(scale_x = c(0, 0), scale_y = c(0, 0), x_name, y_name, extinction_horizon = TRUE, axis_ticks = "black", y_line = "black", rev_ax = TRUE) {
  theme <- theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 23, angle = 90, color = axis_ticks),
                 axis.title.x = element_text(vjust = 0.5, size = 26),
                 axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 23, colour = "black"),
                 axis.title.y = element_text(size = 26),
                 plot.title = element_text(size = 23, face="bold", hjust = 0.5),
                 strip.text.x = element_text(angle = 70), # rotates facet labels
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = y_line),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())
  
  extinction_horizon_gg <- geom_vline(xintercept = 0, colour = "brown1", linetype = "dashed") # horizontal dashed line at extinction horizon
  
  rev_ax_gg <- coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") # reverse axes
  
  parts <- list(theme)
  
  if(rev_ax == TRUE) {
    parts <- c(parts, rev_ax_gg)
  }

  if(class(scale_x) == "numeric") {
    scale_x <- scale_x_continuous(expand = scale_x, name = x_name)
    parts <- c(parts, scale_x)
  }
  
  if(class(scale_y) == "numeric") {
    scale_y <- scale_y_continuous(expand = scale_y, name = y_name)
  } else {
    scale_y <- scale_y_continuous(name = y_name)
  }
  parts <- c(parts, scale_y)
  
  if(extinction_horizon == TRUE) {
    parts <- c(parts, extinction_horizon_gg)
  }
  
  plot_common_parameters <- parts
}

