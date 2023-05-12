library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

#abundance <- read_csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.csv")
meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")

raw_abundance <- read_delim("~/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/phil_running_3/1_species/summarised_filtered_species_otu.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# define a function to calculate the proportion
rel <- function(x) {(x / sum(x)) * 100}
# apply the rel function across each column starting with a 3
relative_abundance <- raw_abundance %>% mutate(across(starts_with('3'), rel))
#relative_abundance <- relative_abundance %>% t() %>% as.data.frame() %>% setNames(relative_abundance[,1])

relative_abundance <- relative_abundance %>% t()  %>% as.data.frame() 
names(relative_abundance) <- as.character(unlist(relative_abundance[1,]))
# remove the first row
relative_abundance <- relative_abundance[-1,]

# add the row names as a column
relative_abundance <- tibble::rownames_to_column(relative_abundance, "ID") %>% as.data.frame()



relative_abundance <- relative_abundance %>% left_join(dplyr::select(meta, c('ID', 'Group', 'Country', 'Age', 'Antibiotics_taken_before_sampling_yes_no_assumptions')), by = 'ID')

abundance <- relative_abundance %>% select('ID', 'Group', 'Country', 'Age', "Antibiotics_taken_before_sampling_yes_no_assumptions", "Prevotella copri", 'Romboutsia timonensis')
abundance$`Prevotella copri` <- as.numeric(as.character(abundance$`Prevotella copri`))
abundance$`Romboutsia timonensis` <- as.numeric(as.character(abundance$`Romboutsia timonensis`))


abundance$age_range <- cut(abundance$Age,
                       breaks=c(0, 5.99, 15.99, Inf),
                       labels=c('0-5', '6-15', '15+'))


abundance <- filter(abundance, Country != 'Nepal')
abundance <- filter(abundance, Group != 'Carrier')
# abundance <- rename(abundance, ProportionPrevotellacopri = `Proportion Prevotella copri`)

abundance_malawi <- filter(abundance, Country == 'Malawi')
abundance_bangladesh <- filter(abundance, Country == 'Bangladesh')

abundance_bangladesh_no_outlier <- filter(abundance_bangladesh, `Prevotella copri` < 10)

give.n <- function(x){
  return(c(y = 17, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}




# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
# https://stackoverflow.com/questions/21468380/overlay-geom-points-on-geom-boxplotfill-group


give.n <- function(x){
  return(c(y = 12, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}

b_pc <- ggplot(abundance_bangladesh_no_outlier, aes(x = age_range, y = `Prevotella copri`, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage P. copri') + 
  ggtitle('Bangladesh') +
  theme(legend.position="none") + 
  stat_compare_means(label = "p.format", label.y = 10) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) 

b_pc

b_pc_sp <- ggplot(abundance_bangladesh, aes(x = Age, y = `Prevotella copri`)) + 
  geom_point(aes(colour = factor(Group))) +
  ggtitle('Bangladesh') +
  theme(legend.position="none")

b_pc_sp

give.n <- function(x){
  return(c(y = 17, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}

b_rt <- ggplot(abundance_bangladesh, aes(x = age_range, y = `Romboutsia timonensis`, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage R. timonensis') + 
  ggtitle('Bangladesh') +
  theme(legend.position="none") + 
  stat_compare_means(label = "p.format", label.y = 15) +
  stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75)) 

b_rt

b_rt_sp <- ggplot(abundance_bangladesh, aes(x = Age, y = `Romboutsia timonensis`)) + 
  geom_point(aes(colour = factor(Group))) +
  ggtitle('Bangladesh') +
  theme(legend.position="none")

b_rt_sp

give.n <- function(x){
  return(c(y = 50, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}


# this is a way to do the comparisons because stat_compare_means balks when some groups have no values.
# do this, then add to the graph as annotate
abundance_malawi_rn <- abundance_malawi %>% rename(prevotella_copri = `Prevotella copri`, romboutsia_timonensis = `Romboutsia timonensis`)
#abundance_malawi_rn <- abundance_malawi %>% rename()

malawi_pc_abundance_compared <- compare_means(prevotella_copri ~ Group, data = abundance_malawi_rn, 
                                           group.by = "age_range")


m_pc <- ggplot(abundance_malawi, aes(x = age_range, y = `Prevotella copri`, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage P. copri') + 
  ggtitle('Malawi') +
  theme(legend.position="none") + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) + 
  annotate("text", x = 3, y = 47, label = paste0('p = ', toString(malawi_pc_abundance_compared[,6]$p.adj)))

m_pc

m_pc_sp <- ggplot(abundance_malawi, aes(x = Age, y = `Prevotella copri`)) + 
  geom_point(aes(colour = factor(Group))) +
  ggtitle('Malawi') +
  theme(legend.position="none")

m_pc_sp


give.n <- function(x){
  return(c(y = 2, label = length(x)))
  #return(data.frame(y = 55, label = paste0('n = ', length(x)))) 
  # experiment with the multiplier to find the perfect position
}

malawi_rt_abundance_compared <- compare_means(romboutsia_timonensis ~ Group, data = abundance_malawi_rn, 
                                              group.by = "age_range")


m_rt <- ggplot(abundance_malawi, aes(x = age_range, y = `Romboutsia timonensis`, fill = factor(Group))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
  ylab('Percentage R. timonensis') + 
  ggtitle('Malawi') +
  theme(legend.position="none") + 
  annotate("text", x = 3, y = 1.7, label = paste0('p = ', toString(malawi_rt_abundance_compared[,6]$p.adj))) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75))

m_rt

m_rt_sp <- ggplot(abundance_malawi, aes(x = Age, y = `Romboutsia timonensis`)) + 
  geom_point(aes(colour = factor(Group)))  +
  ggtitle('Malawi') +
  theme(legend.position="none")

m_rt_sp

(b_pc + m_pc) / (b_rt + m_rt)

(b_pc_sp + m_pc_sp) / (b_rt_sp + m_rt_sp)


# compare the relative abundance of p copri in malawi separated by antibiotic use




# plot a box plot of relative abundance of p copri in malawi separated by antibiotic use

abundance_malawi_cases <- abundance_malawi %>% filter(Group == 'Acute_Typhi')
m_pc_ab <- ggplot(abundance_malawi_cases, aes(x = Antibiotics_taken_before_sampling_yes_no_assumptions, y = `Prevotella copri`)) + geom_boxplot() + stat_compare_means(label = "p.format") + ggtitle('Malawi') + theme(legend.position="none")
m_pc_ab

abundance_malawi_cases_abs <- abundance_malawi_cases %>% filter(Antibiotics_taken_before_sampling_yes_no_assumptions == 'Yes')
abundance_malawi_cases_no_abs <- abundance_malawi_cases %>% filter(Antibiotics_taken_before_sampling_yes_no_assumptions == 'No')
wilcox.test(abundance_malawi_cases_abs$`Prevotella copri`, abundance_malawi_cases_no_abs$`Prevotella copri`)


abundance_bangladesh_cases <- abundance_bangladesh %>% filter(Group == 'Acute_Typhi')
b_pc_ab <- ggplot(abundance_bangladesh_cases, aes(x = Antibiotics_taken_before_sampling_yes_no_assumptions, y = `Prevotella copri`)) + geom_boxplot() + stat_compare_means(label = "p.format") + ggtitle('Bangladesh') + theme(legend.position="none") + ylab('Percentage P. copri')
b_pc_ab

abundance_bangladesh_cases_abs <- abundance_bangladesh_cases %>% filter(Antibiotics_taken_before_sampling_yes_no_assumptions == 'Yes')
abundance_malawi_cases_no_abs <- abundance_bangladesh_cases %>% filter(Antibiotics_taken_before_sampling_yes_no_assumptions == 'No')
wilcox.test(abundance_malawi_cases_abs$`Prevotella copri`, abundance_malawi_cases_no_abs$`Prevotella copri`)



write_csv(abundance, '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/2022.03.28 prevotella copri relative abundance.by_location_and_type.csv')
       