library(readr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

abundance <- read_csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.csv")
meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
abundance <- abundance %>% left_join(select(meta, c('ID', 'Group', 'Country', 'Age')), by = c('Sample' = 'ID'))

abundance$age_range <- cut(abundance$Age,
                       breaks=c(0, 5.99, 15.99, Inf),
                       labels=c('0-5', '6-15', '15+'))


abundance <- filter(abundance, Country != 'Nepal')
abundance <- filter(abundance, Group != 'Carrier')
abundance_malawi <- filter(abundance, Country == 'Malawi')
abundance_bangladesh <- filter(abundance, Country == 'Bangladesh')


ggplot(abundance_malawi, aes(x = Group, y = `Proportion Prevotella copri`)) + geom_boxplot() + ylab('Percentage P. copri') + ggtitle('Malawi')

ggplot(abundance_malawi, aes(x = age_range, y = `Proportion Prevotella copri`, fill = factor(Group))) + geom_boxplot() + ylab('Percentage P. copri') + ggtitle('Malawi')


abundance_bangladesh <- filter(abundance_bangladesh, `Proportion Prevotella copri` < 10)

ggplot(abundance_bangladesh, aes(x = age_range, y = `Proportion Prevotella copri`, fill = factor(Group))) + geom_boxplot() + ylab('Percentage P. copri') + ggtitle('Bangladesh')


ggplot(abundance_bangladesh, aes(x = Group, y = `Proportion Prevotella copri`)) + geom_boxplot() + ylab('Percentage P. copri') + ggtitle('Bangladesh')

write_csv(abundance, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.by_location_and_type.csv')
       