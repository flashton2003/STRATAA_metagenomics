library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

abundance <- read_csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.csv")
meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
abundance <- abundance %>% left_join(dplyr::select(meta, c('ID', 'Group', 'Country', 'Age')), by = c('Sample' = 'ID'))

abundance$age_range <- cut(abundance$Age,
                       breaks=c(0, 5.99, 15.99, Inf),
                       labels=c('0-5', '6-15', '15+'))


abundance <- filter(abundance, Country != 'Nepal')
abundance <- filter(abundance, Group != 'Carrier')
abundance <- rename(abundance, ProportionPrevotellacopri = `Proportion Prevotella copri`)

abundance_malawi <- filter(abundance, Country == 'Malawi')
abundance_bangladesh <- filter(abundance, Country == 'Bangladesh')

abundance_bangladesh_no_outlier <- filter(abundance_bangladesh, ProportionPrevotellacopri < 10)

give.n <- function(x){
  return(c(y = 40, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}




# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# https://stackoverflow.com/questions/21468380/overlay-geom-points-on-geom-boxplotfill-group

b <- ggplot(abundance_bangladesh, aes(x = age_range, y = ProportionPrevotellacopri, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage P. copri') + 
  ggtitle('Bangladesh') +
  theme(legend.position="none") + 
  stat_compare_means(label = "p.format", label.y = 37) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) 

b

#b_no_o <- ggplot(abundance_bangladesh_no_outlier, aes(x = age_range, y = ProportionPrevotellacopri, fill = factor(Group))) + 
#  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
#  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
#  ylab('Percentage P. copri') + 
#  ggtitle('Bangladesh') +
#  theme(legend.position="none") + 
#  stat_compare_means(label = "p.format", label.y = 25) +
#  stat_summary(fun.data = give.n, geom = "text", fun = median,
#               position = position_dodge(width = 0.75))
#b_no_o

malawi_abundance_compared <- compare_means(ProportionPrevotellacopri ~ Group, data = abundance_malawi, 
              group.by = "age_range")


m <- ggplot(abundance_malawi, aes(x = age_range, y = ProportionPrevotellacopri, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage P. copri') + 
  ggtitle('Malawi') +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) + 
  annotate("text", x = 3, y = 52, label = paste0('p = ', toString(malawi_abundance_compared[,6]$p.adj)))

m

b | m

write_csv(abundance, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.by_location_and_type.csv')
       