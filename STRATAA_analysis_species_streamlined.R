library(apeglm)
library(autoplotly)
library(broom)
library(cluster)
library(compositions)
library(dendextend)
library(DESeq2)
library(devtools)
library(dplyr)
library(edgeR)
library(edgeR)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(heatmaply)
library(hrbrthemes)
library(htmlwidgets)
library(MASS)
library(metacoder)
library(microbiome)
library(patchwork)
library(phyloseq)
library(purrr)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape2)
library(rmarkdown)
library(stringr)
library(taxa)
library(tidyverse)
library(vegan)
library(viridis)
library(vsn)
library(zCompositions)
library(VennDiagram)

source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/pairwise_beta.R")
source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/parse_bracken_functions.r")
source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/diff_expr.r")

meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
# i dont think this line does anything, probably just there for historic reasons
names(meta)[names(meta) == "sample_ID"] <- "isolate"
names(meta)[names(meta) == "ID"] <- "isolate"
#meta <- meta %>% mutate(Sex = str_replace(Sex, 'Male  ', 'Male'))

# add age bracket

meta <- meta %>% mutate(age_bracket=cut(Age, breaks=c(0, 1, 5, 15, Inf), labels=c("0-1", "1-5", "6-15", ">15")))


table(meta$Group, meta$Country)
table(meta$Sex, meta$Country)
table(meta$age_bracket, meta$Country)
table(meta$Antibiotics_taken_before_sampling_yes_no_assumptions, meta$Country)

round(prop.table(table(meta$Group, meta$Country)), 2)


## plot of the age per group per country
#ggplot(meta, aes(x = Country, y = Age, fill = Group)) + geom_boxplot()

number_per_country <- meta %>% group_by(Country) %>% summarise(count = n())
number_per_country <- split(number_per_country$count, number_per_country$Country)

eg1 <- meta %>% group_by(Country, Group, Sex) %>% summarise(count = n()) 
for (c in c('Bangladesh', 'Malawi', 'Nepal')) {
  d <- eg1 %>% filter(Country == c)
  p <- ggplot(d, aes(x = Group, y = count, fill = Sex)) + geom_bar(stat ='identity', position = 'fill') + ylab('Proportion') + ggtitle(c)
  show(p)
}

table(meta$Group, meta$Country)

braken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species/"
output_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3/"
taxonomic_level = "species"

run_calc_beta_both_countries <- function(input_braken_folder, out_folder, level){
  species_dir <- file.path(out_folder, '1_species')
  analysis_dir <- file.path(out_folder, '3_analysis')
  beta_dir <- file.path(out_folder, '4_beta')
  if (!dir.exists(out_folder)){ dir.create(out_folder) }
  if (!dir.exists(species_dir)){ dir.create(species_dir) }
  if (!dir.exists(analysis_dir)){ dir.create(analysis_dir) }
  if (!dir.exists(beta_dir)){ dir.create(beta_dir) }
  samples_regexp = "\\d+_\\d+_\\d+" #for patch and strata
  meta <- meta %>% filter(Country == 'Malawi' | Country == 'Bangladesh' | Country == 'Nepal')
  meta$group_country <- paste(meta$Group, meta$Country, sep = '_')
  meta$group_antibiotic <- paste(meta$Group, meta$Antibiotics_taken_before_sampling_yes_no_assumptions, sep = '_')
  
  data_table <- filter_data(input_braken_folder, level, samples_regexp, out_folder)
  
  data_table.unnormalised <- data.matrix(read.csv(paste(out_folder, "/1_species/summarised_filtered_species_otu.txt", sep = ''), header=T, sep = "\t"))
  colnames(data_table.unnormalised) = gsub(pattern = "X", replacement = "", x = colnames(data_table.unnormalised))
  
  #calculate_beta(data_table.unnormalised, meta, quote("country"), out_folder, "country", level, 'all_three')
  calculate_beta(data_table.unnormalised, meta, "Country", out_folder, "Country", level, 'Beta diversity PCoA: species level')
  
}

run_calc_beta_both_countries("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3/", "species")

run_calc_beta_both_countries("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_blantyre_dhaka/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre_dhaka/", "species")

run_calc_beta_one_country <- function(input_braken_folder, out_folder, level, country){
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
  
  calculate_beta(data_table.unnormalised, meta, quote("Country"), out_folder, "Country", level, country)
  group_beta <- calculate_beta(data_table.unnormalised, meta, quote("Group"), out_folder, "Group", level, country)
  sex_beta <- calculate_beta(data_table.unnormalised, meta, quote("Sex"), out_folder, "Sex", level, country)
  age_beta <- calculate_beta(data_table.unnormalised, meta, quote("age_bracket"), out_folder, "age_bracket", level, country)
  
  return(group_beta)
  
}

malawi_group <- run_calc_beta_one_country("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_blantyre/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre/", "species", "Malawi")

nepal_group <- run_calc_beta_one_country("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_kathmandu/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_kathmandu/", "species", "Nepal")

bangladesh_group <- run_calc_beta_one_country("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_dhaka/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_dhaka/", "species", "Bangladesh")


bangladesh_group / malawi_group | nepal_group / plot_spacer()



run_calc_alpha <- function(input_braken_folder, out_folder, level, country){
  print('running calc alpha')
  samples_regexp = "\\d+_\\d+_\\d+" #for patch and strata
  data_table <- filter_data(input_braken_folder, level, samples_regexp, out_folder)
  print('done reading data')
  meta <- meta %>% filter(Country == 'Malawi' | Country == 'Bangladesh')
  View(data_table)
  calculate_alpha(data_table, "kraken_assigned_reads", meta, quote("Group"), out_folder, "Phenotype", level, TRUE)  
}


run_calc_alpha("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre_dhaka/", "species", "Malawi")

#run_calc_alpha("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species/", "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3/", "species", "Malawi")


run_dge <- function(our_metadata, root_dir){
  # filter full_meta
  rownames(our_metadata) <- our_metadata[,1]
  sampledata <- sample_data(our_metadata)
  View(sampledata)
  #OTU
  print(paste(root_dir, "/1_species/summarised_filtered_species_otu.txt", sep = ''))
  otu <- read.csv(paste(root_dir, "/1_species/summarised_filtered_species_otu.txt", sep = ''), header=T, row.names=1, sep = "\t")
  colnames(otu) = gsub(pattern = "X", replacement = "", x = colnames(otu)) 
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  
  #put them in a phyloseq object
  OTU <- t(OTU)
  #print(as.vector(sampledata[,1]))
  #print(length(sampledata))
  #print(length(sampledata[,1]))
  
  #print(sample_names(OTU))
  #print(sample_names(sampledata))
  
  pseq <- phyloseq(OTU, sampledata)
  our_metadata <- meta(pseq)
  otu <- abundances(pseq)
  
  pseq_control_vs_acute <- subset_samples(pseq, Group =="Control_HealthySerosurvey" | Group =="Acute_Typhi")
  subset_meta <- meta(pseq_control_vs_acute)
  View(subset_meta)
  subset_otu <- t(abundances(pseq_control_vs_acute))
  # this is Leo's function for running edgeR GLM.
  result <- glm.edgeR(x=subset_meta$Group, Y=subset_otu, covariates = subset_meta[ , c('Country', 'Sex', 'Age', 'Antibiotics_taken_before_sampling_yes_no_assumptions')])
  #result <- glm.edgeR(x=subset_meta$Group, Y=subset_otu)
  topTags(result, n=10)
  out_folder <- paste(root_dir,'/5_glm/', sep = '')
  if (!dir.exists(out_folder)){ dir.create(out_folder) }
  write.table(topTags(result, n=Inf)$table, file=paste(root_dir,'/5_glm/results_all.country_sex_age_amu.edgeR.tsv', sep = ''),sep='\t',quote=FALSE, col.names=NA)
}


#full_meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, row.names=1, sep = "\t")
run_dge(meta, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3')

bangladesh_meta <- meta %>% filter(Country == 'Bangladesh')
run_dge(bangladesh_meta, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_dhaka')

nepal_meta <- meta %>% filter(Country == 'Nepal')
run_dge(nepal_meta, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_kathmandu')

malawi_meta <- meta %>% filter(Country == 'Malawi')
run_dge(malawi_meta, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre')


combine_and_compare_dges <- function() {
  # dpt is differentially present taxa
  all_sites_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3/5_glm/results_all.country_sex_age_amu.edgeR.tsv'
  bangladesh_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_dhaka/5_glm/results_all.country_sex_age_amu.edgeR.tsv'
  malawi_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre/5_glm/results_all.country_sex_age_amu.edgeR.tsv'
  nepal_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_kathmandu/5_glm/results_all.country_sex_age_amu.edgeR.tsv'
  all_sites_dpt <- read_delim(all_sites_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  all_sites_dpt <- rename(all_sites_dpt, c(species=...1, all_sites_logFC = logFC, all_sites_logCPM = logCPM, all_sites_LR = LR, all_sites_PValue = PValue, all_sites_FDR = FDR))
  
  bangladesh_dpt <- read_delim(bangladesh_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  bangladesh_dpt <- rename(bangladesh_dpt, c(species=...1, bangladesh_logFC = logFC, bangladesh_logCPM = logCPM, bangladesh_LR = LR, bangladesh_PValue = PValue, bangladesh_FDR = FDR))
  bangladesh_dpt_sig <- filter(bangladesh_dpt, bangladesh_FDR <= 0.01)
  bangladesh_dpt_sig_up <- filter(bangladesh_dpt_sig, bangladesh_logFC >= 1)
  bangladesh_dpt_sig_down <- filter(bangladesh_dpt_sig, bangladesh_logFC <= -1)
  malawi_dpt <- read_delim(malawi_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  malawi_dpt <- rename(malawi_dpt, c(species=...1, malawi_logFC = logFC, malawi_logCPM = logCPM, malawi_LR = LR, malawi_PValue = PValue, malawi_FDR = FDR))  
  malawi_dpt_sig <- filter(malawi_dpt, malawi_FDR <= 0.01)
  malawi_dpt_sig_up <- filter(malawi_dpt_sig, malawi_logFC >= 1)
  malawi_dpt_sig_down <- filter(malawi_dpt_sig, malawi_logFC <= -1)
  nepal_dpt <- read_delim(nepal_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  nepal_dpt <- rename(nepal_dpt, c(species=...1, nepal_logFC = logFC, nepal_logCPM = logCPM, nepal_LR = LR, nepal_PValue = PValue, nepal_FDR = FDR))  
  nepal_dpt_sig <- filter(nepal_dpt, nepal_FDR <= 0.01)
  
  combined_dpt <- left_join(all_sites_dpt, bangladesh_dpt, by = "species")
  combined_dpt <- left_join(combined_dpt, malawi_dpt, by = "species")
  combined_dpt <- left_join(combined_dpt, nepal_dpt, by = "species")
  
  write.table(combined_dpt, file='/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.combined_dges.tsv',sep='\t',quote=FALSE, col.names=NA)
  venn.diagram(x = list(bangladesh_dpt$species, malawi_dpt$species, nepal_dpt$species), category.names = c('bangladesh', 'malawi', 'nepal'), filename = '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.venn_diagram.no_filter.png', euler.d = FALSE, scaled = FALSE, height=2200, width=2200)
  venn.diagram(x = list(bangladesh_dpt_sig$species, malawi_dpt_sig$species, nepal_dpt_sig$species), category.names = c('bangladesh', 'malawi', 'nepal'), filename = '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.venn_diagram.fdr_0.01.png', euler.d = FALSE, scaled = FALSE, height=2200, width=2200)
  
  
  venn.diagram(x = list(bangladesh_dpt_sig_up$species, malawi_dpt_sig_up$species), category.names = c('bangladesh', 'malawi'), filename = '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.venn_diagram.fdr_0.01.upreg.png', euler.d = FALSE, scaled = FALSE)
  venn.diagram(x = list(bangladesh_dpt_sig_down$species, malawi_dpt_sig_down$species), category.names = c('bangladesh', 'malawi'), filename = '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.venn_diagram.fdr_0.01.downreg.png', euler.d = FALSE, scaled = FALSE)
  
  
}

combine_and_compare_dges()


