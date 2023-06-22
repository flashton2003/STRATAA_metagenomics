library(dplyr)
#library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)


data_table.raw_counts <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/phil_running_3/1_species/summarised_filtered_species_otu.tsv", header=T, sep = "\t")
names(data_table.raw_counts)[names(data_table.raw_counts) == 'X'] <- 'taxon'
#rownames(data_table.raw_counts) <- data_table.raw_counts$X
#data_table.raw_counts <- subset(data_table.raw_counts, select = -c(X))
colnames(data_table.raw_counts) = gsub(pattern = "X", replacement = "", x = colnames(data_table.raw_counts)) 

data_table <- data_table.raw_counts
y_lab_text <- 'raw counts'


meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
data_table <- pivot_longer(data_table, !taxon)

data_table <- data_table %>% left_join(dplyr::select(meta, c('ID', 'Group', 'Country')), by = c('name' = 'ID'))
healthy_taxa <- c('Prevotella.copri', 'Dialister.sp000434475', 'Prevotella copri', 'Dialister sp000434475')
typhi_taxa <- c( 'Absiella.innocuum', 'Pauljensenia.odontolyticus', 'KLE1796.sp001580115', 'Massiliomicrobiota.timonensis', 'Absiella innocuum', 'Pauljensenia odontolyticus', 'KLE1796 sp001580115', 'Massiliomicrobiota timonensis')
data_table <- filter(data_table, Country != 'Nepal')
data_table <- filter(data_table, Group != 'Carrier')


plot_count <- function(taxon_to_plot, working_data_table){
  to_plot <- filter(working_data_table, taxon == taxon_to_plot)
  View(to_plot)
  p <- ggplot(to_plot, aes(x = taxon, y = value, fill = factor(Group))) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
    ylab(y_lab_text) + 
    xlab('Species')
  p
}


plot_count('Prevotella copri', data_table)
plot_count('Dialister sp000434475', data_table)
plot_count('Absiella innocuum', data_table)
plot_count('Pauljensenia odontolyticus', data_table)
plot_count('KLE1796 sp001580115', data_table)
plot_count('Massiliomicrobiota timonensis', data_table)
