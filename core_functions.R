library(tidyr)
library(readr)
library(magrittr)
library(dplyr)
library(stringr)
library(reshape2)

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
  
  baseline_chars <- meta_subset %>% group_by(Country, Group) %>% summarise(number = n(), median_age = median(Age)) %>% left_join(pct_female, by = c('Group', 'Country')) %>% left_join(pct_antibiotics, by = c('Group', 'Country'))
  return(baseline_chars)
}

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

make_filtered_otu_table <- function(input_braken_folder, taxonomic_level, filename_regex, output_folder, samples_to_include){
  #browser()
  output_folder <- file.path(output_folder, paste('1_', taxonomic_level, sep = ""))
  #create the names of the parsed bracken files - one with all the data and one with the filtered
  summary_file <- file.path(output_folder, paste("summarised_", taxonomic_level, "_kraken.txt", sep = ""))
  otu_file <- file.path(output_folder, paste("summarised_", taxonomic_level, "_otu.txt", sep = ""))
  filtered_otu_file <- file.path(output_folder, paste("summarised_filtered_", taxonomic_level, "_otu.txt", sep = ""))

  #if the filtered data exists, load it and exit
  if (file.exists(filtered_otu_file)){
    filtered_otu <- read.csv(filtered_otu_file, header=T, sep = "\t")
    return(filtered_otu)
  }
  
  #read the bracken output and put it in 1 file
  combined_bracken <- combine_bracken_outputs(braken_folder, filename_regex, samples_to_include, taxonomic_level)
  write.table(combined_bracken, summary_file, row.names = F, quote=FALSE, sep='\t')
  
  #create the otu matrix
  otu <- data.frame(acast(combined_bracken, name ~ sample_ID, value.var = "kraken_assigned_reads", fill=0))
  colnames(otu) = gsub(pattern = "X", replacement = "", x = colnames(otu))
  write.table(otu, otu_file, row.names = T, quote=FALSE, sep='\t')
  
  pdf(file = paste(output_folder, taxonomic_level, "_abundances_hist.pdf", sep = ""))
  hist(as.matrix(otu), 
       breaks = 100, 
       main="Histogram For Unfiltered OTU table", 
       xlab="Relative Abundance", 
       border="blue")
  dev.off()
  
  #old filtering method
  #data <- data %>% filter(fraction_total_reads > 0.002)
  
  #filter: remove OTUs that have non-zero values in <= 10% of samples
  # pa - this is equivalent to removing SNPs with low allele threshold in a GWAS
  filtered_otu <- remove_rare(table=otu, cutoff_pro=0.1)
  colnames(filtered_otu) = gsub(pattern = "X", replacement = "", x = colnames(filtered_otu))
  write.table(filtered_otu, filtered_otu_file, row.names = T, quote=FALSE, sep='\t')
  
  otu_table_rare_removed_norm_cpm <- sweep(filtered_otu, 2, colSums(filtered_otu) , '/')*100
  
  pdf(file = paste(output_folder, taxonomic_level, "_filtered_abundances_hist.pdf", sep = ""))
  hist(as.matrix(otu_table_rare_removed_norm_cpm), 
       main="Histogram For filtered OTU table", 
       breaks = 100, 
       xlab="Relative Abundance", 
       border="blue")
  dev.off()
  
  
  return(filtered_otu)
}


combine_bracken_outputs <- function(braken_folder, filename_regex, samples_to_include, tax_level){
  # compile all reports
  bracken_reports_files <- list.files(braken_folder, full.names = T)
  
  combined_bracken <- NULL
  for (sample in samples_to_include){
    file_handle <- file.path(braken_folder, tax_level, paste(sample, '_R1.bracken', sep = ''))
    print(file_handle)
    file_info <- file.info(file_handle)
    # if the file contains data then add it to the combined table
    if (file_info$size > 99){
      bracken_reports_file <- read_tsv(file_handle)
      #browser()
      #here change the reg expr to match the sample names
      name <- str_extract(rownames(file_info), filename_regex)
      bracken_reports_file$sample_ID <- name
      combined_bracken <- rbind(combined_bracken, bracken_reports_file)
    }
  }
  return(combined_bracken)
}
