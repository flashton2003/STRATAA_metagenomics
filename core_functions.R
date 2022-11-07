
library(readr)
library(magrittr)
library(dplyr)

read_metadata <- function(path_to_metadata){
  meta <- read.csv(path_to_metadata, header=T, sep = "\t")
  # i dont think this line does anything, probably just there for historic reasons
  names(meta)[names(meta) == "sample_ID"] <- "isolate"
  names(meta)[names(meta) == "ID"] <- "isolate"
  meta <- meta %>% mutate(age_bracket=cut(Age, breaks=c(0, 1, 5, 15, Inf), labels=c("0-1", "1-5", "6-15", ">15")))
  
  return(meta)
}

get_baseline_characteristics <- function(meta){
  meta_subset <- meta %>% select(Group, Sex, Country, Age, Antibiotics_taken_before_sampling_yes_no_assumptions)
  
  pct_female <- meta_subset %>% group_by(Group, Country, Sex) %>% summarise(n = n()) %>% pivot_wider(names_from = c(Sex), values_from = n) %>% mutate(pct_fem = (Female / sum(c(Female, Male))) * 100) %>% select(c(Group, Country, pct_fem))
  
  pct_antibiotics <- meta_subset %>% group_by(Group, Country, Antibiotics_taken_before_sampling_yes_no_assumptions) %>% summarise(n = n()) %>% pivot_wider(names_from = c(Antibiotics_taken_before_sampling_yes_no_assumptions), values_from = n)
  pct_antibiotics[is.na(pct_antibiotics)] <- 0
  pct_antibiotics <- pct_antibiotics %>% mutate(pct_anti = (Yes / sum(c(No, Yes, Unknown))) * 100) %>% select(c(Group, Country, pct_anti))
  
  
  baseline_chars <- meta_subset %>% group_by(Group, Country) %>% summarise(median_age = median(Age)) %>% left_join(pct_female, by = c('Group', 'Country')) %>% left_join(pct_antibiotics, by = c('Group', 'Country'))
  return(baseline_chars)
}
