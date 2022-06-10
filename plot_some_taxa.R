library(dplyr)
#library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

'''
1. read in the normalised data and metadata
2. pull out the taxa and samples i want
3. reshape it to be ggplot friendly
4. plot it
'''

data_table.normalised.clr <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/phil_running_3/1_normalized_data/clr_normalised_data.txt", header=T, sep = "\t")
#data_table.normalised.clr <- read_delim("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/phil_running_3/1_normalized_data/clr_normalised_data.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(data_table.normalised.clr) = gsub(pattern = "X", replacement = "", x = colnames(data_table.normalised.clr)) 
data_table.normalised.clr <- as.data.frame(t(data_table.normalised.clr))
colnames(data_table.normalised.clr) = gsub(pattern = "X", replacement = "", x = colnames(data_table.normalised.clr)) 
data_table.normalised.clr$taxon <- rownames(data_table.normalised.clr)
rownames(data_table.normalised.clr) <- NULL




data_table <- data_table.normalised.clr
y_lab_text <- 'CLM transformed reads'


meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
data_table <- pivot_longer(data_table, !taxon)

data_table <- data_table %>% left_join(dplyr::select(meta, c('ID', 'Group', 'Country')), by = c('name' = 'ID'))
healthy_taxa <- c('Prevotella.copri', 'Dialister.sp000434475', 'Prevotella copri', 'Dialister sp000434475')
typhi_taxa <- c( 'Absiella.innocuum', 'Pauljensenia.odontolyticus', 'KLE1796.sp001580115', 'Massiliomicrobiota.timonensis', 'Absiella innocuum', 'Pauljensenia odontolyticus', 'KLE1796 sp001580115', 'Massiliomicrobiota timonensis')
data_table <- filter(data_table, Country != 'Nepal')
data_table <- filter(data_table, Group != 'Carrier')
healthy_all <- filter(data_table, taxon %in% healthy_taxa)
typhi_all <- filter(data_table, taxon %in% typhi_taxa)

p_healthy_all <- ggplot(healthy_all, aes(x = taxon, y = value, fill = factor(Group))) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab(y_lab_text) + xlab('Species')

p_typhi_all <- ggplot(typhi_all, aes(x = taxon, y = value, fill = Group)) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab(y_lab_text) + xlab('Species')

p_healthy_all / p_typhi_all + plot_layout(heights = c(2, 4))


malawi_healthy <- filter(healthy_all, Country == 'Malawi')
malawi_typhi <- filter(typhi_all, Country == 'Malawi')

p_healthy_malawi <- ggplot(malawi_healthy, aes(x = taxon, y = value, fill = factor(Group))) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab('CLM transformed reads') + xlab('Species')
p_typhi_malawi <- ggplot(malawi_typhi, aes(x = taxon, y = value, fill = Group)) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab('CLM transformed reads') + xlab('Species')

p_healthy_malawi / p_typhi_malawi + plot_layout(heights = c(2, 4))



bangladesh_healthy <- filter(healthy_all, Country == 'Bangladesh')
bangladesh_typhi <- filter(typhi_all, Country == 'Bangladesh')

p_healthy_bangladesh <- ggplot(bangladesh_healthy, aes(x = taxon, y = value, fill = factor(Group))) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab('CLM transformed reads') + xlab('Species')
p_typhi_bangladesh <- ggplot(bangladesh_typhi, aes(x = taxon, y = value, fill = Group)) + geom_boxplot() + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + coord_flip() + ylab('CLM transformed reads') + xlab('Species')

p_healthy_bangladesh / p_typhi_bangladesh + plot_layout(heights = c(2, 4))


plot_count <- function(taxon_to_plot, working_data_table){
  View(working_data_table)
  to_plot <- filter(working_data_table, taxon == taxon_to_plot)
  View(to_plot)
  p <- ggplot(to_plot, aes(x = taxon, y = value, fill = factor(Group))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + ylab(y_lab_text) + xlab('Species')
  p
}


plot_count('Prevotella copri', data_table)
plot_count('Dialister sp000434475', data_table)
plot_count('Absiella innocuum', data_table)
plot_count('Pauljensenia odontolyticus', data_table)
plot_count('KLE1796 sp001580115', data_table)
plot_count('Massiliomicrobiota timonensis', data_table)


