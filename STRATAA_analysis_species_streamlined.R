library(broom)
library(tidyverse)
library(stringr)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(heatmaply)
library(ggplot2)
library(readr)
library(ggfortify)
library(cluster)
library(autoplotly)
library(htmlwidgets)
library(compositions)
library(zCompositions)
library(vegan)
library(dendextend)
library(dplyr)
library(ggpubr)
library(edgeR)
library(rmarkdown)
library(gridExtra)
library(microbiome)
library(phyloseq)
library(apeglm)
library(vsn)
library(DESeq2)
library(edgeR)
library(metacoder)
library(taxa)
library(purrr)
library(ggrepel)
library(devtools)
library(hrbrthemes)
library(MASS)


source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/bin/parse_bracken_functions.r")
source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/bin/pairwise_beta.R")
source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/bin/diff_expr.r")


meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/0_metadata/full_meta.txt", header=T, sep = "\t")
# i dont think this line does anything, probably just there for historic reasons
names(meta)[names(meta) == "sample_ID"] <- "isolate"
names(meta)[names(meta) == "ID"] <- "isolate"


braken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species/"
output_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/phil_running_3/"
taxonomic_level = "species"

run_calc_beta <- function(input_braken_folder, out_folder, level, country){
  species_dir <- file.path(out_folder, '1_species')
  analysis_dir <- file.path(out_folder, '3_analysis')
  beta_dir <- file.path(out_folder, '4_beta')
  
  if (!dir.exists(out_folder)){ dir.create(out_folder) }
  if (!dir.exists(species_dir)){ dir.create(species_dir) }
  if (!dir.exists(analysis_dir)){ dir.create(analysis_dir) }
  if (!dir.exists(beta_dir)){ dir.create(beta_dir) }
  
  #give filename patern as a reg exp
  samples_regexp = "\\d+_\\d+_\\d+" #for patch and strata
  meta <- meta %>% filter(Country == country)
  
  data_table <- filter_data(input_braken_folder, level, samples_regexp, out_folder)
  
  data_table.unnormalised <- data.matrix(read.csv(paste(out_folder, "/1_species/summarised_filtered_species_otu.txt", sep = ''), header=T, sep = "\t"))
  colnames(data_table.unnormalised) = gsub(pattern = "X", replacement = "", x = colnames(data_table.unnormalised))
  
  calculate_beta(data_table.unnormalised, meta, quote("Country"), out_folder, "Country", level)
  calculate_beta(data_table.unnormalised, meta, quote("Group"), out_folder, "Group", level)
}

run_calc_beta("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_blantyre/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/phil_blantyre/", "species", "Malawi")

run_calc_beta("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_kathmandu/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome_from_Leo/Leonardos_analysis/phil_kathmandu/", "species", "Nepal")

