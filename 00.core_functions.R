library(dplyr)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(microbiome)
library(phyloseq)
library(readr)
library(reshape2)
library(rlist)
library(stringr)
library(tidyr)
library(vegan)
library(VennDiagram)
library(Maaslin2)
library(forestplot)
library(RColorBrewer)
library(kableExtra)
library(patchwork)
library(forcats)

#source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/config.R")

read_metadata <- function(path_to_metadata){
  meta <- read.csv(path_to_metadata, header=T, sep = "\t")
  # View(meta)
  # i dont think this line does anything, probably just there for historic reasons
  names(meta)[names(meta) == "sample_ID"] <- "isolate"
  names(meta)[names(meta) == "ID"] <- "isolate"
  meta <- meta %>% rename(SampleID = Lane)
  meta <- meta %>% mutate(age_bracket=cut(Age, breaks=c(0, 1, 5, 15, Inf), labels=c("0-1", "1-5", "6-15", ">15")))
  meta$group_country <- paste(meta$Group, meta$Country, sep = '_')
  meta$group_antibiotic <- paste(meta$Group, meta$Antibiotics_taken_before_sampling_assumptions, sep = '_')
  
  # keep only one sample per participant, the one with the most reads.
  meta <- meta %>% arrange(desc(number_of_reads)) %>% distinct(StudyID, .keep_all = TRUE)
  meta <- meta %>% filter(!is.na(isolate) & isolate != "")
  # remove the carriers from bangladesh
  meta <- meta %>% filter(!(Group == 'Carrier' & Country == 'Bangladesh'))
  return(meta)
}


get_baseline_characteristics <- function(meta){
  meta_subset <- meta %>% select(Group, Sex, Country, Age, Antibiotics_taken_before_sampling_assumptions)
  pct_female <- meta_subset %>% group_by(Group, Country, Sex) %>% summarise(n = n()) %>% pivot_wider(names_from = c(Sex), values_from = n) %>% mutate(pct_fem = (Female / sum(c(Female, Male))) * 100) %>% select(c(Group, Country, pct_fem))
  pct_antibiotics <- meta_subset %>% group_by(Group, Country, Antibiotics_taken_before_sampling_assumptions) %>% summarise(n = n()) %>% pivot_wider(names_from = c(Antibiotics_taken_before_sampling_assumptions), values_from = n)
  pct_antibiotics[is.na(pct_antibiotics)] <- 0
  pct_antibiotics <- pct_antibiotics %>% mutate(pct_anti = (Yes / sum(c(No, Yes))) * 100) %>% select(c(Group, Country, pct_anti))
  
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


run_make_filtered_otu_table <- function(output_folder, tax_level, countries, filename_regex, braken_folder, meta){
  species_dir <- file.path(output_folder, '1_species')
  alpha_dir <- file.path(output_folder, '3_alpha')
  
  if (!dir.exists(output_folder)){ dir.create(output_folder) }
  if (!dir.exists(species_dir)){ dir.create(species_dir) }
  if (!dir.exists(alpha_dir)){ dir.create(alpha_dir) }
  
  samples_to_include <- meta %>% filter(Country %in% countries) %>% select('isolate')
  samples_to_include <- samples_to_include$isolate
  #View(samples_to_include)
  filtered_otu.unnormalised <- make_filtered_otu_table(braken_folder, tax_level, filename_regex, output_folder, samples_to_include)
  # not needed at the moment, but leaving in just in case
  #colnames(filtered_otu.unnormalised) = gsub(pattern = "X", replacement = "", x = colnames(data_table.unnormalised))
  return(filtered_otu.unnormalised)
}


make_filtered_otu_table <- function(braken_folder, taxonomic_level, filename_regex, output_folder, samples_to_include){
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
  
  pdf(file = file.path(output_folder, paste(taxonomic_level, "_abundances_hist.pdf", sep = "")))
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
  
  pdf(file = file.path(output_folder, paste(taxonomic_level, "_filtered_abundances_hist.pdf", sep = ""))) 
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


calculate_beta <- function(data, meta, output_folder){
  #make the first column headers
  rownames(meta) <- meta[,1]
  #meta <- meta %>% remove_rownames %>% column_to_rownames()
  
  #and order
  meta <- meta[ order(row.names(meta)), ]
  # need to transpose the OTU table from species as rows to species as columns
  data <- as.matrix(t(data))
  #order data table as well to be sure
  data <- data[ order(row.names(data)), ]
  rownames(data) <- gsub("#","_",rownames(data))
  rownames(data) <- gsub("X","",rownames(data))
  #make sure that you have the same ids
  meta_names <- rownames(meta)
  data_names <- rownames(data)
  common.ids <- intersect(rownames(meta), rownames(data))
  #browser()
  #make sure you have the correct meta
  meta <- meta[common.ids,]
  data <- data[common.ids,]
  
  #transform the table to relative abundances
  data <- sweep(data, 1, rowSums(data),'/')
  #calculate bray curtis distance matrix
  d.bray <- vegdist(data)
  
  #transform it to a matrix to save it
  beta_matrix <- as.matrix(d.bray)
  #View(beta_matrix)
  #Perform PCoA
  pc.bray <- cmdscale(d.bray, k=4, eig = T)
  pcoa.var <- round(pc.bray$eig/sum(pc.bray$eig)*100, 1)
  pcoa.values <- pc.bray$points
  pcoa.data <- data.frame(Sample = rownames(pcoa.values), X=pcoa.values[,1], Y = pcoa.values[,2])
  #pcoa.data <- data.frame(X=pcoa.values[,1], Y = pcoa.values[,2])
  
  if (!dir.exists(output_folder)){ dir.create(output_folder) }
  
  #save the table
  pairwise_out_file <-  file.path(output_folder, "pairwise_beta.txt")
  first_2d_coords_out_file <- file.path(output_folder, "first_2d_coords.txt")
  pcoa_var_file <- file.path(output_folder, "pcoa_var.txt")
  write_delim(as.data.frame(beta_matrix), pairwise_out_file)
  #write.table(beta_matrix, pairwise_out_file, row.names=T, col.names=T, sep = "\t")
  write_delim(as.data.frame(pcoa.data), first_2d_coords_out_file)
  #write.table(pcoa.data, first_2d_coords_out_file, col.names = T, sep = "\t")
  write_delim(as.data.frame(pcoa.var), pcoa_var_file)
  output <- list(pcoa.data = pcoa.data, pcoa.var = pcoa.var)
  saveRDS(d.bray, file.path(output_folder, 'd.bray.RDS'))
  saveRDS(meta, file.path(output_folder, 'meta_for_permanova.RDS'))
  
}


plot_beta <- function(pcoa.data, pcoa.var, to_plot){
  #plot and save
  #file_path <- paste(output_folder, prefix, "_beta_PCoA.pdf", sep = "")
  #get the metadata column to paint the plot
  output_plots <- list()
  if ("Country" %in% to_plot) {
    country_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = Country)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Country")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    
    
    output_plots <- list.append(output_plots, "country_plot" = country_plot)
  }
  
  if ("Group" %in% to_plot) {
    group_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = Group)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Group")) +
      geom_point() #+
      #coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    output_plots <- list.append(output_plots, "group_plot" = group_plot)
  }
  
  if ("Sex" %in% to_plot) {
    sex_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = Sex)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Sex")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    output_plots <- list.append(output_plots, "sex_plot" = sex_plot)
  }
  
  if ("age_bracket" %in% to_plot) {
    age_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = age_bracket)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Age bracket")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    output_plots <- list.append(output_plots, "age_plot" = age_plot)
  }
  
  if ("group_country" %in% to_plot) {
    group_country_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = group_country)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Group/Country")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    output_plots <- list.append(output_plots, "group_country_plot" = group_country_plot)
  }
  
  if ("group_antibiotic" %in% to_plot) {
    group_antibiotic_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = group_antibiotic)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="Group/Antibiotic")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
    output_plots <- list.append(output_plots, "group_antibiotic_plot" = group_antibiotic_plot)
  }
    
  if ("sequencing_lane" %in% to_plot) {
    sequencing_lane_plot <- ggplot(pcoa.data, aes(x=X, y=Y, colour = sequencing_lane)) + 
      xlab(paste("MDS1 - ", pcoa.var[[1]][1], "%", sep="")) + 
      ylab(paste("MDS2 - ", pcoa.var[[1]][2], "%", sep="")) + 
      guides(colour=guide_legend(title="sequencing_lane")) +
      geom_point() +
      coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  output_plots <- list.append(output_plots, "sequencing_lane" = sequencing_lane_plot)
  }
  
  
  #+
    #theme(legend.position="none")
  #+ geom_text(aes(label=Sample),hjust=0, vjust=0)
  #g1 <- ggplot(pcoa.data, aes(x=X, y=Y)) + 
  #  ggtitle(title) + xlab(paste("MDS1 - ", pcoa.var[1], "%", sep="")) + ylab(paste("MDS2 - ", pcoa.var[2], "%", sep="")) + 
  #  geom_point() #+ geom_text(aes(label=Sample),hjust=0, vjust=0)
  #ggsave(file_path)
  #print(g1)
  return(output_plots)
}

run_beta_diversity <- function(metaphlan_data, metadata, groups_of_interest){
  metaphlan_data_mat <- metaphlan_data |> 
    as.matrix() |>
    t()

  dist_mat <- vegdist(metaphlan_data_mat)
  cmd_res <- cmdscale(dist_mat, 
                    k = (nrow(metaphlan_data_mat) - 1),
                    eig = TRUE)
  pcoa_df <- tibble(PC1 = cmd_res$points[,1], 
                  PC2 = cmd_res$points[,2],
                  PC3 = cmd_res$points[,3],
                  PC4 = cmd_res$points[,4])

  pcoa_df <- mutate(pcoa_df, SampleID = rownames(metaphlan_data_mat)) %>% left_join(metadata, by = 'SampleID')
  # pn <- adonis(dist_mat~Sex*Group*Age*Antibiotics_taken_before_sampling_assumptions, data = metadata, permutations = 100)
  # View(pn)
  if ('High Vi-titre' %in% groups_of_interest){
    pn <- adonis2(dist_mat~Sex*Group*Age, data = metadata, permutations = 1000)
  } else if ('patch_baseline' %in% groups_of_interest){
    pn <- adonis2(dist_mat~Gender*Diagnosis*age_at_challenge, data = metadata, permutations = 1000)
  }
  else {
    pn <- adonis2(dist_mat~Sex*Group*Age*Antibiotics_taken_before_sampling_assumptions, data = metadata, permutations = 1000)
  }
  # View(pn)
  # pn_with_var_names <- cbind(rownames(pn$aov.tab), data.frame(pn$aov.tab, row.names = NULL))
  # pn_res <- pn_with_var_names %>% rename(variable = `rownames(pn$aov.tab)`) %>% mutate(is_it_significant = ifelse(`Pr..F.` < 0.01, 'significant', 'not_significant')) %>% arrange(desc(is_it_significant), desc(R2)) %>% select(!c(Df, SumsOfSqs, MeanSqs, F.Model))
  pn_res <- pn %>% mutate(is_it_significant = ifelse(`Pr(>F)` < 0.01, 'significant', 'not_significant')) %>% arrange(`Pr(>F)`, desc(is_it_significant), desc(R2)) %>% select(!c(Df, SumOfSqs, F))
  # View(pn_res)
  p <- list(pcoa_df = pcoa_df, pn_res = pn_res)
  return(p)
}

strataa_metaphlan_beta <- function(metaphlan_data, metadata, countries_of_interest, groups_of_interest, participant_group_colours){
  # View(metadata)
  # View(metaphlan_data)
  metadata_coi <- metadata %>% filter(Country %in% countries_of_interest) %>% filter(Group %in% groups_of_interest)
  metadata_coi_ids <- metadata_coi %>% pull(SampleID) 
  metaphlan_data_coi <- metaphlan_data %>% select(all_of(metadata_coi_ids))
  # View(metadata_coi)
  # View(metaphlan_data_coi)
  rbd_output <- run_beta_diversity(metaphlan_data_coi, metadata_coi, groups_of_interest)
  # View(rbd_output)
  # View(rbd_output$pcoa_df)
  # View(participant_group_colours)
  title <- paste(countries_of_interest, collapse = "_")
  pc12 <- ggplot(rbd_output$pcoa_df, aes(x = PC1, y = PC2, colour = Group)) + 
    geom_point() +
    ggtitle(title) +
    # coord_fixed() + 
    scale_color_manual(values = participant_group_colours) +
    stat_ellipse()
  show(pc12)
  pc34 <- ggplot(rbd_output$pcoa_df, aes(x = PC3, y = PC4, colour = Group)) + 
    geom_point() +
    ggtitle(title) +
    # coord_fixed() + 
    scale_color_manual(values = participant_group_colours) +
    stat_ellipse()
  p <- list(pc12 = pc12, pc34 = pc34, pn_res = rbd_output$pn_res)
  return(p)
}


metaphlan_alpha <- function(metaphlan_data, metaphlan_metadata, countries_of_interest, groups_of_interest, comparisons, participant_group_colours) {
  metaphlan_data <- metaphlan_data %>% select(!lowest_taxonomic_level)
  # View(metaphlan_data)
  alpha <- rbiom::alpha.div(as.matrix(metaphlan_data))
  alpha_meta <- left_join(alpha, metaphlan_metadata, by = c('Sample' = 'SampleID')) %>% filter(Country %in% countries_of_interest) %>% filter(Group %in% groups_of_interest)
  # View(alpha_meta)
  # if the length of countries of interest is greater than 1, then include Country in the model
  # if carrier is in the groups of interest, don't include antibiotics in the model (because carriers and healthy are all assumed "No" for antibiotics)
  if (length(countries_of_interest) > 1) {
    if ('High Vi-titre' %in% groups_of_interest){
      alpha_anova <- aov(Shannon ~ Country * Sex *Group * Age, data = alpha_meta)
    } else {
      alpha_anova <- aov(Shannon ~ Country * Sex *Group * Age * Antibiotics_taken_before_sampling_assumptions, data = alpha_meta)
    }
    
  } else if (length(countries_of_interest) == 1) {
    if ('High Vi-titre' %in% groups_of_interest){
      alpha_anova <- aov(Shannon ~ Sex *Group * Age, data = alpha_meta)
    } else {
      alpha_anova <- aov(Shannon ~ Sex *Group * Age * Antibiotics_taken_before_sampling_assumptions, data = alpha_meta)
    }
  }
  
  alpha_anova_summary <- summary(alpha_anova)
  
  alpha_anova_summary_with_var_names <- cbind(rownames(alpha_anova_summary[[1]]), data.frame(alpha_anova_summary[[1]], row.names = NULL))
  alpha_anova_summary_with_var_names <- alpha_anova_summary_with_var_names %>% mutate(is_it_significant = ifelse(`Pr..F.` < 0.01, 'significant', 'not_significant')) %>% arrange(desc(is_it_significant), `Pr..F.`)

  # this makes a list of vectors, each of which contains a pair of countries to compare
  if (length(countries_of_interest) > 1){
    country_comps <- combn(countries_of_interest, 2)
    country_comps <- split(country_comps, col(country_comps))
    alpha_by_country <- ggplot(alpha_meta, aes(x=Country, y=Shannon)) + 
      labs(x="", y="Shannon diversity") +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold')) +
      stat_compare_means(size = 4, label = "p.format", comparisons = country_comps, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf), symbols = c("***", "**", "*", "ns")))

  } else {
    alpha_by_country <- NULL
  }
  # View(alpha_meta)
  alpha_meta$Group <- factor(alpha_meta$Group, levels = c("Household contact", "Acute typhoid", "High Vi-titre"))

  alpha_plot_group <- ggboxplot(alpha_meta, facet.by = "Country", y = "Shannon", x = "Group", color = "Group") + 
    stat_compare_means(comparisons = comparisons, label = 'p.signif', symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf), symbols = c("***", "**", "*", "ns"))) + 
    rremove("x.text") + 
    rremove("xlab") + 
    rremove("x.ticks") +
    scale_color_manual(values = participant_group_colours)
  # +rotate_x_text(angle = 45)
  show(alpha_plot_group)
  
  antibiotic_comps <- list(c('Yes', 'No'))

  alpha_plot_antibiotics <- ggboxplot(alpha_meta, facet.by = "Country", y = "Shannon", x = "Antibiotics_taken_before_sampling_assumptions", color = "Antibiotics_taken_before_sampling_assumptions") + 
    stat_compare_means(comparisons = antibiotic_comps, label = 'p.signif', symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf), symbols = c("***", "**", "*", "ns"))) + 
    rremove("x.text") + 
    rremove("xlab") + 
    rremove("x.ticks") # +rotate_x_text(angle = 45)

    # stat_compare_means(comparisons = antibiotic_comps, label = 'p.signif', symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, Inf), symbols = c("***", "**", "*", "ns"))) + 

  df <- data.frame(country = c(rep("Bangladesh", 6), rep("Malawi", 6), rep("Nepal", 6)),
                  antibiotic_presentation = rep(c("Pre", "Post"), 9),
                  age_bracket = rep(c("0-4", "5-14", "15-24"), each = 6),
                  value = rnorm(18))
  
if ('Acute typhoid' %in% groups_of_interest){
  # Create a box plot for each combination of country, antibiotic_presentation, and age_bracket
  
  alpha_meta_no_unknowns <- alpha_meta %>% filter(Antibiotics_taken_before_sampling_assumptions != 'Unknown')
  
  p <- ggplot(alpha_meta_no_unknowns, aes(x = age_bracket , y = Shannon, fill = Antibiotics_taken_before_sampling_assumptions)) +
    geom_boxplot() +
    facet_grid(cols = vars(Country )) +
    labs(title = "Box plot by country, antibiotic_presentation, and age_bracket")
  show(p)
}

results <- list(alpha_by_country = alpha_by_country, alpha_anova_summary_with_var_names = alpha_anova_summary_with_var_names, alpha_plot_group = alpha_plot_group, alpha_plot_antibiotics = alpha_plot_antibiotics)
return(results)

}


calculate_alpha <- function(data, meta, group, output_folder, prefix, inc_country){
  
  #calculate alpha diverities 
  #input_matrix <- acast(data, sample_ID ~ name, value.var = summary_column, fill=0)
  #View(data)
  data <- t(data)
  
  
  #make the first column headers
  rownames(meta) <- meta[,1]
  #meta <- meta %>% remove_rownames %>% column_to_rownames()
  
  #and order
  meta <- meta[ order(row.names(meta)), ]
  
  #order data table as well to be sure
  data <- data[ order(row.names(data)), ]
  rownames(data) <- gsub("#","_",rownames(data))
  rownames(data) <- gsub("X","",rownames(data))
  #make sure that you have the same ids
  meta_names <- rownames(meta)
  data_names <- rownames(data)
  #View(meta_names)
  #View(data_names)
  common.ids <- intersect(rownames(meta), rownames(data))
  #View(common.ids)
  #browser()
  #make sure you have the correct meta
  meta <- meta[common.ids,]
  data <- data[common.ids,]
  
  alpha <- vegan::diversity(data, index="shannon")
  
  if (isTRUE(inc_country)) {
    s <- summary(aov(alpha ~ Country * Sex *Group * age_bracket * Antibiotics_taken_before_sampling_assumptions, data = meta))
  } else {
    s <- summary(aov(alpha ~ Sex *Group * age_bracket * Antibiotics_taken_before_sampling_assumptions, data = meta))
  }
  
  

  
  alpha_table <- tibble(isolate = names(alpha), alpha = alpha)
  #View(alpha_table)
  #alpha_table <- alpha_table %>% mutate(isolate = substring(isolate, 2))
  #View(meta)
  #View(alpha_table)
  #alpha_table <- left_join(meta, alpha_table, by="isolate")
  #View(alpha_table)
  
  #for t.test
  #all pairwise combinations
  #print(group)
  #print(meta[,eval(group)])
  #print(levels(meta[,eval(group)]))
  #my_comparisons <- combn(levels(meta[,eval(meta$Group)]), 2, simplify = F)
  
  #pairwise t test
  f <- paste("alpha~", group, sep = "")
  #pv <- compare_means(as.formula(f),  data = alpha_table, method = "t.test")
  
  #gr <- pv$p <= 0.05
  
  write_delim(alpha_table, file.path(output_folder, "3_alpha", "shannon.alpha_results.tsv"))
  #output_folder <- paste(output_folder, "3_alpha/", sep = "")
  #if (!dir.exists(output_folder)){ dir.create(output_folder) }
  
  title <- paste("Alpha div.(T test): ", prefix, sep = "")
  
  
  alpha_table <- alpha_table %>% filter(!is.na(alpha))

  
  return(s)
}


calculate_dge <- function(our_metadata, output_folder, tax_level, otu, covariates, groups_for_comparison){
  # manipulate full_meta
  rownames(our_metadata) <- our_metadata[,1]
  sampledata <- sample_data(our_metadata)

  #otu <- read.csv(otu_file, header=T, row.names=1, sep = "\t")
  colnames(otu) = gsub(pattern = "X", replacement = "", x = colnames(otu)) 
  OTU = otu_table(otu, taxa_are_rows = TRUE)
  
  #put them in a phyloseq object
  OTU <- t(OTU)
  pseq <- phyloseq(OTU, sampledata)
  our_metadata <- meta(pseq)
  otu <- abundances(pseq)
  
  # need to filter out the samples that aren't included in the comparison
  # this is the old way of doing it, but subset_samples can't take variables as arguments (!)
  #pseq_control_vs_acute <- subset_samples(pseq, Group == group1_for_comp | Group == group2)
  # see below link for old bug report and work around
  # https://github.com/joey711/phyloseq/issues/335
  # take the pseq sample data
  var_values <- sample_data(pseq)
  # make a boolean of the length sample_data(pseq) that evaluates whether the Group is in the groups to be compared
  var_bool <- var_values$Group %in% groups_for_comparison
  # use this boolean to prune the pseq object samples
  # i've checked, and this gives equivalent results to the previous way of doing it.
  pseq_control_vs_acute <- prune_samples(var_bool, pseq)
  subset_meta <- meta(pseq_control_vs_acute)
  #View(subset_meta)
  subset_otu <- t(abundances(pseq_control_vs_acute))
  # this is Leo's function for running edgeR GLM.
  result <- glm_edgeR(x=subset_meta$Group, Y=subset_otu, covariates = subset_meta[ , covariates])
  #result <- glm.edgeR(x=subset_meta$Group, Y=subset_otu)
  topTags(result, n=10)
  dge_out_folder <- file.path(output_folder,'5_glm')
  if (!dir.exists(dge_out_folder)){ dir.create(dge_out_folder) }
  covar_initials <- paste(str_sub(covars, 1, 6), sep = '', collapse = '')
  write.table(topTags(result, n=Inf)$table, file=file.path(dge_out_folder, paste('results_all', paste(groups_for_comparison[1], 'vs', groups_for_comparison[2], sep = '_'), covar_initials, 'edgeR.tsv', sep = '.')), sep='\t',quote=FALSE, col.names=NA)
}


glm_edgeR <- function(x, Y, covariates=NULL,use.fdr=TRUE, estimate.trended.disp=TRUE,verbose=TRUE) {
  # x is the independent variable
  # Y is a matrix of samples x dependent variables
  # x is the Group
  # y is the OTU
  # returns p-values
  # drop samples that are NA for, by default, the Group (i.e. acute typhi/healthy control/etc)
  ix <- !is.na(x)
  Y <- Y[ix,]
  x <- x[ix]
  # if covariates have been given
  if(!is.null(covariates)){
    # if the dimensions of them are null then it isn't a matrix/data frame, so convert it to one
    if(is.null(dim(covariates))){
      covariates <- as.data.frame(covariates)
    }
    # drop samples that are NA for the independent variable
    covariates <- covariates[ix,,drop=F]
    
    # drop constant covariates
    covariates <- covariates[,apply(covariates,2,function(xx) length(unique(xx)) > 1),drop=F]
  }
  
  if(verbose) cat('Making DGEList...\n')
  d <- DGEList(count=t(Y), group=x)
  View(d)
  if(verbose) cat('calcNormFactors...\n')
  d <- calcNormFactors(d, method = "TMM" )
  if(!is.null(covariates)){
    covariates <- as.data.frame(covariates)
    # combines the Group (i.e. acute typhi etc) back with the co-variates
    covariates <- cbind(x, covariates)
    covariates <- droplevels(covariates)
    design <- model.matrix(~ ., data=covariates)
    # this makes the below
    # x <- data.frame(c('case', 'control', 'carrier'), c(5, 17, 8))
    # model.matrix(~ ., data = x)
    #  (Intercept) c..case....control....carrier..case c..case....control....carrier..control c.5..17..8.
    # 1           1                                   1                                      0           5
    # 2           1                                   0                                      1          17
    # 3           1                                   0                                      0           8
  } else {
    design <- model.matrix(~x)
  }
  
  if(verbose) cat('estimate common dispersion...\n')
  d <- estimateGLMCommonDisp(d, design)
  if(estimate.trended.disp){
    if(verbose) cat('estimate trended dispersion...\n')
    d <- estimateGLMTrendedDisp(d, design)
  }
  if(verbose) cat('estimate tagwise dispersion...\n')
  d <- estimateGLMTagwiseDisp(d,design)
  
  if(verbose) cat('fit glm...\n')
  fit <- glmFit(d,design)
  if(verbose) cat('likelihood ratio test...\n')
  lrt <- glmLRT(fit,coef=2)
  
  return(lrt)
}


combine_and_compare_edgeRs <- function(to_combine, location_names){
  #View(to_combine)
  print(length(to_combine))
  if (length(to_combine) == 2) {
    combined_dpt <- left_join(to_combine[[1]], to_combine[[2]], by = 'species', suffix = location_names[1:2])
  }
  else {
    print('combine_and_compare_edgeRs only setup for combining 2 dfs right now')
    quit()
  }
  #if (length(to_combine) == 3) {
  #  combined_dpt <- left_join(to_combine[[1]], to_combine[[2]], by = 'species', suffix = location_names[1:2])
  #  combined_dpt <- left_join(combined_dpt, to_combine[[3]], by = 'species', suffix = c(NA, location_names[3]))
  #}
  
  # across selects all the columns that start with FDR, and combined with the filter, selects only rows with FDR < 0.05
  # i'm not really sure what the ~ .x is about though?
  sig <- combined_dpt %>% filter(across(starts_with('FDR'), ~ .x < 0.05))
  sig_up <- sig %>% filter(across(starts_with('logFC'), ~ .x >= 1))
  sig_down <- sig %>% filter(across(starts_with('logFC'), ~ .x <= -1))
  View(sig_up)
  View(sig_down)
  #venn.diagram(x = list(bangladesh_dpt$species, malawi_dpt$species, nepal_dpt$species), category.names = c('Bangladesh', 'Malawi', 'Nepal'), filename = file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'venn_diagram.no_filter.png', sep = '.')), euler.d = FALSE, scaled = FALSE, height=2200, width=2200)
  #print(paste('Jaccard of unfiltered taxa of Bang, Mal, Nep = ', jaccard(bangladesh_dpt$species, malawi_dpt$species, nepal_dpt$species)), sep = '')
  output <- list(sig_up = sig_up, sig_down = sig_down)
  return(output)
}


make_name <- function(output_folder, covars, comp){
  dge_out_folder <- file.path(output_folder,'5_glm')
  covar_initials <- paste(str_sub(covars, 1, 6), sep = '', collapse = '')
  name=file.path(dge_out_folder, paste('results_all', comp, covar_initials, 'edgeR.tsv', sep = '.'))
  return(name)
}


jaccard <- function(input_1, input_2, input_3){
  intersection = length(intersect(intersect(input_1, input_2), input_3))
  union = length(union(union(input_1, input_2), input_3))
  return (intersection/union)
}

# add a fucntion make_clean to take in the maaslin_prevalence and 
# make it into a clean table

make_clean <- function(maaslin_prevalence, variable) {
  # Filter the maaslin_prevalence by variable
  maaslin_filtered <- maaslin_prevalence %>% 
    filter(variable == {{ variable }})
  
  # Select the columns that aren't NA for that variable
  maaslin_clean <- maaslin_filtered %>% 
    select_if(~ !all(is.na(.)))
  # View(maaslin_clean)
  # Combine the last four columns into a single string column
  maaslin_clean <- maaslin_clean %>% 
    unite("prevalence", -c(1:2), sep = ";", na.rm = TRUE)
  # View(maaslin_clean)
  # Return the clean data frame
  return(maaslin_clean)
}


run_make_clean <- function(maaslin_prevalence) {
  # this function runs make_clean for different variables
  # and then combines the output into a single data frame
  group_clean <- make_clean(maaslin_prevalence, variable = "Group")
  sex_clean <- make_clean(maaslin_prevalence, variable = "Sex")
  antibiotics_clean <- make_clean(maaslin_prevalence, variable = "Antibiotics_taken_before_sampling_assumptions")
  clean_output <- rbind(group_clean, sex_clean, antibiotics_clean)
  View(clean_output)
  return(clean_output)
}


add_prevalence_to_maaslin_output <- function(output_dir, prevalence) {
  # Load the MaAsLin output file
  maaslin_output <- read_tsv(file.path(output_dir, "all_results.tsv"))
  # rename the feature values in the maaslin output and the prevalence (from the metaphlan output) - annoying, not sure why this needed to be done.
  maaslin_output <- maaslin_output %>% mutate(feature = str_replace_all(feature, "\\.", "_"))
  prevalence <- prevalence %>% mutate(feature = str_replace_all(feature, "\\|", "_"))
  
  # Pivot the prevalence data frame to make it wider - each row now contains the median prevalence for each feature in each sample for each variable
  prevalence_wider <- prevalence %>% 
    pivot_wider(names_from = value, values_from = Median_Prevalence)
  # View(prevalence_wider)
  # get the clean prevalence data frame, which contains the prevalence in the two different groups for each variable separated by semi-colons
  clean_prevalence <- run_make_clean(prevalence_wider)
  # Join the MaAsLin output file with the metadata
  maaslin_prevalence <- maaslin_output %>% 
    left_join(clean_prevalence, by = join_by("feature" == "feature", "metadata" == "variable"))
  # View(maaslin_prevalence)
  
  # Write the prevalence data frame to a file
  write_tsv(maaslin_prevalence, file.path(output_dir, "all_results.prevalence.tsv"))
}


calculate_prevalence <- function(feature_data, metadata, country, groups_for_analysis, strataa_patch) {
  # Filter the metadata to get the correct country and groups for analysis
  if (strataa_patch == 'patch'){
    metadata_filtered <- metadata
  } else if (strataa_patch == 'strataa') {
    metadata_filtered <- metadata %>% 
    filter(Country == country, Group %in% groups_for_analysis)
  }
  # set SampleID column instead of rownames.
  metadata_filtered <- metadata_filtered %>% 
    mutate(SampleID = rownames(metadata_filtered))
  rownames(metadata_filtered) <- NULL
  # for the feature data as well
  feature_data <- feature_data %>% 
    mutate(feature = rownames(feature_data))
  rownames(feature_data) <- NULL
  # feature data is in a square matrix, so pivot it to long format
  feature_data <- feature_data %>% 
    pivot_longer(cols = -feature, names_to = "SampleID", values_to = "prevalence")
  # View(feature_data)
  # View(metadata_filtered)
  feature_data <- feature_data %>% mutate(SampleID = str_replace(SampleID, '#', '_'))
  metadata_filtered <- metadata_filtered %>% mutate(SampleID = str_replace(SampleID, '#', '_'))
  # Join the feature data with the metadata
  # remove everything with NA for group, these are things from the feature data that aren't in the filtered metadata
  # View(feature_data)
  # View(metadata_filtered)
  feature_metadata <- feature_data %>% 
    left_join(metadata_filtered, by = "SampleID") %>% filter(!is.na(Group))
  # View(feature_metadata)
  # Remove the columns that aren't needed
  feature_metadata <- feature_metadata %>% mutate(sequencing_lane = NULL, Age = NULL, Country = NULL)
  # pivot feature_metadata to long format, before, was multiple metadata columns per row, after this will be one metadata column per row
  feature_metadata_longer <- feature_metadata %>% 
    pivot_longer(cols = -c(feature, SampleID, prevalence), names_to = "variable", values_to = "value")
  
  # View(feature_metadata_longer)
  # variable is e.g. Group/Sex, value is e.g. AcuteTyphoid/Male etc.
  prevalence <- feature_metadata_longer %>% 
    group_by(feature, variable, value) %>% 
    summarize(Median_Prevalence = median(prevalence), )
  
  # View(prevalence)
  # Return the prevalence data frame
  return(prevalence)
}

run_maaslin <- function(feature_data, metadata, output_root, country, groups_for_analysis, variables_for_analysis, norm, trans, reference_groups, input_type){
  ifelse(!dir.exists(output_root), dir.create(output_root), FALSE)
  metadata_to_analyse <- metadata %>% filter(Country == country, Group %in% groups_for_analysis) %>% filter(Antibiotics_taken_before_sampling_assumptions %in% c('Yes', 'No'))
  # View(metadata_to_analyse)
  # View(unique(metadata_to_analyse$Group))
  # View(unique(metadata_to_analyse$Sex))
  # View(dim(feature_data))
  # View(dim(feature_data %>% na.omit()))
  # View(ncol(feature_data))
  # prevalence is a data frame where each row is the prevalence of a feature in a sample, and there is a row for each metadata group that a sample belongs to
  # e.g. feature 1, group 1, prevalence 0.5
  #      feature 1, group 2, prevalence 0.5
  
  vars_for_dirname <- paste(variables_for_analysis, collapse = '.')
  output_dir <- file.path(output_root, paste(country, paste(groups_for_analysis, collapse = '_vs_'), vars_for_dirname, sep = '_'))
  Maaslin2(input_data = feature_data, input_metadata = metadata_to_analyse, analysis_method = "LM", min_prevalence = 0,
           normalization  = norm,
           transform = trans,
           output         = output_dir, 
           fixed_effects  = variables_for_analysis,
           reference = reference_groups)
  if (input_type == 'metaphlan'){
    prevalence <- calculate_prevalence(feature_data, metadata, country, groups_for_analysis, 'strataa')
    add_prevalence_to_maaslin_output(output_dir, prevalence)
  }
  
}


read_in_maaslin <- function(country, groups_to_analyse, variables_for_analysis, type_of_input){
  if (type_of_input == 'metaphlan'){
    root_folder <- maaslin_taxonomic_output_root_folder
  } else if (type_of_input == 'bigmap'){
    root_folder <- maaslin_functional_output_root_folder
  }
  vars_for_dirname <- paste(variables_for_analysis, collapse = '.')
  maaslin_output_dir <- file.path(root_folder, paste(country, paste(groups_to_analyse, collapse = '_vs_'), vars_for_dirname, sep = '_'))
  maaslin_results <- read_delim(file.path(maaslin_output_dir, "all_results.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  if (type_of_input == 'metaphlan'){
      maaslin_results$lowest_taxonomic_level <- sapply(str_split(maaslin_results$feature, "\\."), function(x) x[length(x)])
  maaslin_results <- maaslin_results %>% relocate(lowest_taxonomic_level, .after = feature)
  } else if (type_of_input == 'bigmap'){
    maaslin_results <- maaslin_results %>% separate_wider_delim(feature, delim = 'Entryname.', names_sep = '', cols_remove = FALSE) %>% separate_wider_delim(feature2, delim = '..OS.', names_sep = '') %>% separate_wider_delim(feature22, delim = '..SMASH', names_sep = '') %>% dplyr::select(!c(feature1, feature222)) %>% rename(MGC_class = feature21, Species = feature221, feature = featurefeature)
  }
  return(maaslin_results)
}


filter_taxonomic_maaslin <- function(maaslin_results){
  maaslin_results_species <- maaslin_results %>% filter(grepl('s__', feature)) %>% filter(!grepl('t__', feature))
  maaslin_results_species_group <- maaslin_results_species %>% filter(metadata == "Group")
  return(maaslin_results_species_group)
}


basic_maaslin_stats <- function(taxonomic_maaslin_filtered, country, variables_for_analysis, vars_for_dirname){
  # View(maaslin_results_species_group)
  vars_for_dirname <- paste(variables_for_analysis, collapse = '.')
  volcano_plot <- ggplot(aes(x = coef, y = -log10(qval)), data = taxonomic_maaslin_filtered) + 
    geom_point() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggtitle(paste(country, variables_for_analysis, 'covars=', vars_for_dirname, sep = '_'))
  # show(volcano_plot)
  maaslin_results_sig <- taxonomic_maaslin_filtered %>% filter(qval < 0.05)
  basic_maaslin_results <- list(volcano_plot = volcano_plot, maaslin_results_sig = maaslin_results_sig)
  return(basic_maaslin_results)
  # geom_text(aes(label = feature), size = 2, vjust = 1, hjust = 1) + 
  # View(maaslin_results_clean_sig)
  # sig_per_metadata <- maaslin_results_clean_sig %>% group_by(metadata) %>% summarize(n_sig = n())
  # View(sig_per_metadata)
}


inner_join_maaslins <- function(first_set_maaslin_results, second_set_maaslin_results, first_suffix, second_suffix, type_of_input){
  # thanks chatgpt!
  if (type_of_input == 'metaphlan'){
    # we include the lowest taxonomic level in the join so that it doesn't get duplicated as bang/mal
    combined_df <- first_set_maaslin_results %>%
      inner_join(second_set_maaslin_results, by = c("feature", 'metadata', 'value', 'lowest_taxonomic_level'), suffix = c(first_suffix, second_suffix))
  } else if (type_of_input == 'bigmap'){
    combined_df <- first_set_maaslin_results %>%
      inner_join(second_set_maaslin_results, by = c("MGC_class", "Species", "feature", 'metadata', 'value'), suffix = c(first_suffix, second_suffix))
  }
  return(combined_df)
}


summarise_mgcs <- function(functional_maaslin_results){
  func_maaslin_filtered <- functional_maaslin_results %>% filter(metadata == 'Group', value == 'Control_HealthySerosurvey', qval < 0.05)
  func_maaslin_filtered <- func_maaslin_filtered %>% mutate(direction_of_assc = ifelse(coef > 0, 'assc_health', 'assc_disease'))
  func_maaslin_summarised <- func_maaslin_filtered %>% group_by(MGC_class, direction_of_assc) %>% summarise(n = n())
  return(func_maaslin_summarised)
}

filter_combined_maaslins <- function(combined_df){
  # Filter the combined data frame based on the conditions for coef > 0
  # View(combined_df)
  filtered_df_positive_coef <- combined_df %>%
    filter(if_all(starts_with('qval'), ~ .x < 0.05)) %>% filter(if_all(starts_with('coef'), ~ .x > 0))
  
  # Filter the combined data frame based on the conditions for coef < 0
  filtered_df_negative_coef <- combined_df %>%
    filter(if_all(starts_with('qval'), ~ .x < 0.05)) %>% filter(if_all(starts_with('coef'), ~ .x < 0))
  
  # Return the two filtered data frames
  return(list(positive_coef = filtered_df_positive_coef, negative_coef = filtered_df_negative_coef, all_features = combined_df))
}


run_inner_join_maaslins <- function(maaslin_results, suffixes, variables_for_output_name, groups_to_analyse, analysis_type, maaslin_output_root_folder){
  # maaslin_results is a list of dataframes
  # suffixes is a list of suffixes for each dataframe, to be applied to the output of the join
  groups_for_dirname <- paste(groups_to_analyse, collapse = '_vs_')
  vars_for_output_dirname <- paste(variables_for_output_name, collapse = '.')
  
  if (length(maaslin_results) >= 2){
    combined_maaslins <- inner_join_maaslins(maaslin_results[[1]], maaslin_results[[2]], suffixes[[1]], suffixes[[2]], analysis_type)
  }

  if (length(maaslin_results) == 3){
    combined_maaslins <- inner_join_maaslins(combined_maaslins, maaslin_results[[3]], '', suffixes[[3]], analysis_type)
    # rename coef, std_err, N, N.not.0, pval, qval columns to include the last suffix
    combined_maaslins <- combined_maaslins %>% rename(!!paste('coef', suffixes[[3]], sep = '') := coef, !!paste('stderr', suffixes[[3]], sep = '') := stderr, !!paste('N', suffixes[[3]], sep = '') := N, !!paste('N.not.0', suffixes[[3]], sep = '') := N.not.0, !!paste('pval', suffixes[[3]], sep = '') := pval, !!paste('qval', suffixes[[3]], sep = '') := qval)
    # combined_maaslins 
  }

  if (length(maaslin_results) >3) {
    print('run_inner_join_maaslins only setup for combining 3 dfs right now')
    quit()
  }
  
  
  
  # bang_maaslin_only <- bang_maaslin %>% filter(qval < 0.05) %>% filter(!feature %in% combined_maaslins_positive_coef$feature) %>% filter(!feature %in% combined_maaslins_negative_coef$feature) %>% arrange(desc(coef))
  # mwi_maaslin_only <- malawi_maaslin %>% filter(qval < 0.05) %>% filter(!feature %in% combined_maaslins_positive_coef$feature) %>% filter(!feature %in% combined_maaslins_negative_coef$feature) %>% arrange(desc(coef))
  # View(bang_maaslin_only)
  # print(nrow(bang_maaslin_only) + nrow(combined_maaslins_positive_coef) + nrow(combined_maaslins_negative_coef))
  # print(nrow(bang_maaslin %>% filter(qval < 0.05)))
  
  # combined_maaslins_positive_coef %>%  kbl() %>% kable_styling()
  # write_csv(combined_maaslins_positive_coef, file.path(maaslin_output_root_folder, paste(groups_for_dirname, vars_for_output_dirname, 'combined_maaslins_positive_coef.csv', sep = '.')))
  # combined_maaslins_negative_coef %>%  kbl() %>% kable_styling()
  # write_csv(combined_maaslins_negative_coef, file.path(maaslin_output_root_folder, paste(groups_for_dirname, vars_for_output_dirname, 'combined_maaslins_negative_coef.csv', sep = '.')))
  # write_csv(bang_maaslin_only, file.path(maaslin_output_root_folder, paste(groups_for_dirname, vars_for_output_dirname, 'bangladesh_only.csv', sep = '.')))
  # write_csv(mwi_maaslin_only, file.path(maaslin_output_root_folder, paste(groups_for_dirname, vars_for_output_dirname, 'malawi_only.csv', sep = '.')))

  # combined_results <- list(positive_coef = combined_maaslins_positive_coef, negative_coef = combined_maaslins_negative_coef, bang_maaslin_only = bang_maaslin_only, mwi_maaslin_only = mwi_maaslin_only)
  # combined_results <- list(positive_coef = combined_maaslins_positive_coef, negative_coef = combined_maaslins_negative_coef)
  return(combined_maaslins)
}


run_combine_edgeR <- function(groups_for_comparison, bangladesh_covars, malawi_covars){
  comp <- paste(groups_for_comparison[1], 'vs', groups_for_comparison[2], sep = '_')
  bang_covar_initials <- paste(str_sub(bangladesh_covars, 1, 6), sep = '', collapse = '')
  mal_covar_initials <- paste(str_sub(malawi_covars, 1, 6), sep = '', collapse = '')
  
  blantyre_output_path <- file.path(output_folder_blantyre, '5_glm', paste('results_all', comp, mal_covar_initials, 'edgeR.tsv', sep = '.'))
  dhaka_output_path <- file.path(output_folder_dhaka, '5_glm', paste('results_all', comp, bang_covar_initials, 'edgeR.tsv', sep = '.'))
  
  blantyre_output <- read_delim(blantyre_output_path, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  dhaka_output <- read_delim(dhaka_output_path, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  
  # this is the ridiculous shit you have to do in R to get a list of dataframes.
  to_combine_paths <- c(blantyre_output_path, dhaka_output_path)
  # this is a vec of strings that will be passed to dge_output to give sensible names to the joined output
  # needs to be same order as to_combine_paths
  to_combine_strings_for_join <- c("_mal", "_bang")
  to_combine <- list()
  for (i in seq_along(to_combine_paths)) {
    to_combine[[i]] <- read_delim(to_combine_paths[i], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
    to_combine[[i]] <- to_combine[[i]] %>% rename(species=...1)
  }
  
  dge_output <- combine_and_compare_edgeRs(to_combine, to_combine_strings_for_join)
  
  #options(scipen = 999)
  #options(scipen = 0) # to switch scientific notation back on
  
  sig_up_for_writing <- dge_output$sig_up %>% select(species, logFC_bang, FDR_bang, logFC_mal, FDR_mal) 
  sig_up_for_writing %>% kbl() %>% kable_styling()
  write_csv(sig_up_for_writing, file.path(combined_output_root, 'dge', paste(the_date, comp, covar_initials, 'sig_up_dge.csv', sep = '.')))
  
  sig_down_for_writing <- dge_output$sig_down %>% select(species, logFC_bang, FDR_bang, logFC_mal, FDR_mal)
  sig_down_for_writing %>% kbl() %>% kable_styling()
  write_csv(sig_down_for_writing, file.path(combined_output_root, 'dge', paste(the_date, comp, covar_initials, 'sig_down_dge.csv', sep = '.')))
  #covar_initials <- paste(str_sub(covars, 1, 6), sep = '', collapse = '')
  #combined_dge_output_folder <- file.path(combined_output_root, paste(covar_initials, 'combined_edgeR'))
}


plot_species_of_interest <- function(prevalence_meta, species_of_interest, country_of_interest, groups_of_interest, participant_group_colours){
  
  # we specify the annotation positions manually using these list of lists
  country_species_y_value <- list('Bangladesh'=list('s__Prevotella_copri_clade_A'=6, 's__Clostridium_SGB6179'=5, 's__GGB4266_SGB5809' = 14, 's__Haemophilus_parainfluenzae'=1.5, 's__Romboutsia_timonensis'=7), 'Malawi'=list('s__Prevotella_copri_clade_A'=47, 's__Clostridium_SGB6179'=0.6, 's__GGB4266_SGB5809'=5.5, 's__Haemophilus_parainfluenzae'=3, 's__Romboutsia_timonensis'=2, 's__Lachnospiraceae_bacterium'=2), 'Nepal'=list('s__Prevotella_copri_clade_A'=20, 's__Clostridium_SGB6179'=1, 's__GGB4266_SGB5809'=5, 's__Haemophilus_parainfluenzae'=0.25, 's__Romboutsia_timonensis'=2))
  # for the number of samples in each group, we need to use this function, and i can't figure out how to dynamically return 
  # the y value dependiing on the input from within the stat_summary call so we have to do this shit show.
  # country_species_fun_data <- list('Bangladesh'=list('s__Prevotella_copri_clade_A'=function(x){return(c(y = 6, label = NA))}, 's__Clostridium_SGB6179'=function(x){return(c(y = 4.5, label = length(x)))}, 's__GGB4266_SGB5809' = function(x){return(c(y = 10, label = length(x)))}, 's__Haemophilus_parainfluenzae'=function(x){return(c(y = 1.25, label = length(x)))}, 's__Romboutsia_timonensis'=function(x){return(c(y = 6, label = length(x)))}, 's__Lachnospiraceae_bacterium'=function(x){return(c(y = 1.75, label = length(x)))}), 'Malawi'=list('s__Prevotella_copri_clade_A'=function(x){return(c(y = 45, label = NA))}, 's__Clostridium_SGB6179'=function(x){return(c(y = 0.5, label = length(x)))}, 's__GGB4266_SGB5809'=function(x){return(c(y = 5, label = length(x)))}, 's__Haemophilus_parainfluenzae'=function(x){return(c(y = 2.5, label = length(x)))}, 's__Romboutsia_timonensis'=function(x){return(c(y = 1.5, label = length(x)))}, 's__Lachnospiraceae_bacterium'=function(x){return(c(y = 1.5, label = length(x)))}), 'Nepal'=list('s__Prevotella_copri_clade_A'=function(x){return(c(y = 6, label = NA))}, 's__Clostridium_SGB6179'=function(x){return(c(y = 4.5, label = length(x)))}, 's__GGB4266_SGB5809' = function(x){return(c(y = 10, label = length(x)))}, 's__Haemophilus_parainfluenzae'=function(x){return(c(y = 1.25, label = length(x)))}, 's__Romboutsia_timonensis'=function(x){return(c(y = 6, label = length(x)))}, 's__Lachnospiraceae_bacterium'=function(x){return(c(y = 1.75, label = length(x)))}))
  # country_species_fun_data <- list()
  # give.n <- country_species_fun_data[[country_of_interest]][[species_of_interest]]
  
  prevalence_meta_filtered <- prevalence_meta %>% filter(lowest_taxonomic_level == species_of_interest) %>% filter(Country == country_of_interest) %>% filter(Group %in% groups_of_interest)  

  abundance_compared <- compare_means(prevalence ~ Group, data = prevalence_meta_filtered, group.by = "age_bracket") %>% arrange(age_bracket)

  p <- ggplot(prevalence_meta_filtered, aes(x = age_bracket, y = prevalence, fill = factor(Group))) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
    ylab(paste('Percentage', species_of_interest, sep = ' ')) + 
    ggtitle(country_of_interest) +
    theme(legend.position="none") + 
    scale_color_manual(values = participant_group_colours) +
    scale_fill_manual(values = participant_group_colours) #+
    # stat_summary(fun.data = country_species_fun_data[[country_of_interest]][[species_of_interest]], geom = "text", fun = median, position = position_dodge(width = 0.75))
    # stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75))

  # add the p values to the plot, has to be done like this because there are
  # not enough comparisons for most of the groups, so do manually and add
  # in. there are more groups with enough to compare in bangladesh, hence
  # the multiple annotations.
  if (country_of_interest %in% c('Malawi', 'Nepal')){
    p <- p + annotate("text", x = 3, y = country_species_y_value[[country_of_interest]][[species_of_interest]], label = paste0('p = ', toString(abundance_compared[,6]$p.adj)))
  } else if (country_of_interest == 'Bangladesh' & species_of_interest != 's__Prevotella_copri_clade_A'){
    p <- p + 
      annotate("text", x = 3, y = country_species_y_value[[country_of_interest]][[species_of_interest]], label = paste0('p = ', toString(abundance_compared[,6]$p.adj[[2]]))) +
      annotate("text", x = 4, y = country_species_y_value[[country_of_interest]][[species_of_interest]], label = paste0('p = ', toString(abundance_compared[,6]$p.adj[[3]])))
  }
    
  # because there is an outlier for s__Lachnospiraceae_bacterium in Malawi, we use ggforce facet zoom to show a zoomed in version alongside the full version.
  if (species_of_interest == 's__Lachnospiraceae_bacterium' & country_of_interest == 'Malawi'){
    p <- p + ggforce::facet_zoom(ylim = c(0, 3))
    return(p)
  } 
  else if (species_of_interest == 's__Prevotella_copri_clade_A' & country_of_interest == 'Bangladesh') {
    # p <- p + ggforce::facet_zoom(ylim = c(0, 6.5))
    p <- ggplot(prevalence_meta_filtered, aes(x = age_bracket, y = prevalence, fill = factor(Group))) + 
      geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
      geom_point(position = position_jitterdodge(), aes(group=factor(Group), colour = factor(Group))) +
      ylab(paste('Percentage', species_of_interest, sep = ' ')) + 
      ggtitle(country_of_interest) +
      theme(legend.position="none") + 
      scale_color_manual(values = participant_group_colours) +
      scale_fill_manual(values = participant_group_colours) +
      ylim(c(0, 10)) #+
      # stat_summary(fun.data = country_species_fun_data[[country_of_interest]][[species_of_interest]], geom = "text", fun = median, position = position_dodge(width = 0.75)) +
   p <- p + 
      annotate("text", x = 3, y = country_species_y_value[[country_of_interest]][[species_of_interest]], label = paste0('p = ', toString(abundance_compared[,6]$p.adj[[2]]))) +
      annotate("text", x = 4, y = country_species_y_value[[country_of_interest]][[species_of_interest]], label = paste0('p = ', toString(abundance_compared[,6]$p.adj[[3]])))
    return(p)
  } 
    else if (species_of_interest == 's__GGB4266_SGB5809' & country_of_interest == 'Bangladesh') {
    p <- p + ggforce::facet_zoom(ylim = c(0, 17))
    return(p)
  } else {
    return(p)
  }
}


run_plot_species_of_interest <- function(prevalence_meta, species_of_interest, participant_group_colours){
  # we do this so that we can combine plots for the same species from bang and mal together
  m <- plot_species_of_interest(strataa_metaphlan_data_longer_meta, species_of_interest, 'Malawi', c('Acute typhoid', 'Household contact'), participant_group_colours)
  m <- m + theme(text = element_text(size = 15))
  b <- plot_species_of_interest(strataa_metaphlan_data_longer_meta, species_of_interest, 'Bangladesh', c('Acute typhoid', 'Household contact'), participant_group_colours)
  # show(m)
  b <- b + theme(text = element_text(size = 15))
  n <- plot_species_of_interest(strataa_metaphlan_data_longer_meta, species_of_interest, 'Nepal', c('Acute typhoid', 'Household contact'), participant_group_colours)
  # show(m)
  n <- n + theme(text = element_text(size = 15))
  
  # to combine the y-axis titles, need to make another plot which is just text (eyeroll)
  p4 <- ggplot(data.frame(l = m$labels$y, x = 1, y = 1)) +
      geom_text(aes(x, y, label = l), angle = 90, size = 6) + 
      theme_void() +
      coord_cartesian(clip = "off")
  # set all the y axis titles as blank
  b$labels$y <- m$labels$y <- n$labels$y <- " "
  # combine the y  axis title plot with the three plots
  p <- p4 + (m / b / n) + plot_layout(widths = c(1, 15))#& theme(

  outhandle_name <- paste(species_of_interest, 'Acute typhoid',  'Household contact', 'Bangladesh', 'Malawi', 'Nepal', 'png', sep = '.')
  ggsave(file.path(maaslin_taxonomic_output_root_folder, outhandle_name), p, width = 10, height = 7)
  show(p)
  return(p)
}


plot_per_country_abundance <- function(phyla_clean_metadata, country, group_order){
    phyla_clean_country <- phyla_clean_metadata %>% filter(Country == country) %>% filter(Group %in% group_order)
    # View(phyla_clean_country)
    phyla_clean_country_fct <- phyla_clean_country %>%
      mutate(Group = factor(Group, levels = group_order),  # Convert group to a factor with the desired order
      group_order_numeric = as.numeric(Group),  # Create a new numeric variable based on the order of group
      sample = fct_reorder(sample, group_order_numeric))  # Reorder sample based on group_order_numeric
    # View(phyla_clean_country_fct)
    p <- ggplot(data = phyla_clean_country_fct, aes(x = sample, y = relative_abundance, fill = clade_name)) + 
        geom_bar(stat = "identity") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) + 
        scale_x_discrete(breaks=phyla_clean_country_fct$sample, labels=phyla_clean_country_fct$Group) + 
        ggtitle(country) +
        ylim(0, 100.01)
    # show(p)
    return(p)
}


prep_data_to_plot_phyla <- function(strataa_metaphlan_data, metadata_select){
  # get the taxa that are phyla, not classes, or below (species etc), and tidy the data.
  phyla <- strataa_metaphlan_data %>% mutate(clade_name = rownames(strataa_metaphlan_data)) %>% filter(grepl("p__", clade_name)) %>% filter(!grepl("c__", clade_name)) %>% pivot_longer(!c(clade_name, lowest_taxonomic_level), names_to = "sample", values_to = "relative_abundance")

  # relative_abundance > 1 returns a list of TRUE/FALSE values, which is then summed to get the number of samples in which the phylum is present at > 1% relative abundance.
  # then we filter to only keep phyla that are present at 1% in at least 10% of samples.
  phyla_to_exclude <- phyla %>% group_by(clade_name) %>% 
      summarise(count = sum(relative_abundance > 1)) %>% 
      filter(count < (length(unique(phyla$sample)) / 10)) %>% 
      pull(clade_name)
  # View(phyla_to_exclude)

  # in order to make each sample add up to 100, we need to add the excluded taxa back in as a single "rare taxa" phylum.
  # first we need to calculate the relative abundance of the excluded taxa in each sample.
  excluded_phyla <- phyla %>%
    filter(clade_name %in% phyla_to_exclude) %>% group_by(sample) %>% summarise(relative_abundance = sum(relative_abundance))
  # then make a column with the name "rare_taxa" for each sample, annd bind it to the excluded taxa data.
  rare_taxa_column <- data.frame(lowest_taxonomic_level = c(rep("rare_taxa", nrow(excluded_phyla))), clade_name = c(rep("rare_taxa", nrow(excluded_phyla))))
  excluded_phyla <- cbind(rare_taxa_column, excluded_phyla)

  # then remove the excluded taxa from the phyla data.
  phyla_clean <- phyla %>%
    filter(!(clade_name %in% phyla_to_exclude))
  # and add the excluded taxa back in.
  phyla_clean <- rbind(phyla_clean, excluded_phyla)
  # View(phyla_clean)
  # View(excluded_phyla)

  # colnames(strataa_metaphlan_metadata)
  # metadata_select <- strataa_metaphlan_metadata %>% dplyr::select(SampleID, Group, Country)
  phyla_clean_metadata <- phyla_clean %>% left_join(metadata_select, by = c("sample" = "SampleID"))
  # phyla_clean_metadata <- phyla_clean_metadata %>% mutate(Group = ifelse(Group == "Acute typhoid", "Typhi", Group)) %>% mutate(Group = ifelse(Group == "Control", "Healthy", Group))
  return(phyla_clean_metadata)
}


get_maaslin_results_for_species <- function(combined_maaslins, species_of_interest){
  # grepl for the species name, and then filter for the Group metadata (which is the case/control assc) and then select all columns except feature, metadata, value
  # maaslin_results_for_species <- combined_maaslins %>% filter(grepl(species_of_interest, lowest_taxonomic_level)) %>% filter(metadata == 'Group') %>% select(!c(feature, metadata, value))
  maaslin_results_for_species <- combined_maaslins %>% filter(lowest_taxonomic_level %in% species_of_interest) %>% filter(metadata == 'Group') %>% select(!c(feature, metadata, value))
  # pivot the data so that the columns are now the country, and the values are the coef, stderr, N, N.not.0, pval, qval
  maaslin_results_for_species <- maaslin_results_for_species %>%
    pivot_longer(
      cols = starts_with("coef") | starts_with("stderr") | starts_with("N") | starts_with("N.not.0") | starts_with("pval") | starts_with("qval"),
      names_to = c(".value", "country"),
      names_pattern = "(.*)_(.*)"
    )
  # change the country names to the full names
  maaslin_results_for_species$country <- ifelse(maaslin_results_for_species$country == 'bang', 'Bangladesh', maaslin_results_for_species$country)
  maaslin_results_for_species$country <- ifelse(maaslin_results_for_species$country == 'mal', 'Malawi', maaslin_results_for_species$country)
  maaslin_results_for_species$country <- ifelse(maaslin_results_for_species$country == 'nep', 'Nepal', maaslin_results_for_species$country)
  maaslin_results_for_species$country <- ifelse(maaslin_results_for_species$country == 'patch', 'CHIM', maaslin_results_for_species$country)
  # change the country to a factor, and set the levels to the desired order
  desired_order <- c("Bangladesh", "Malawi", "Nepal", "CHIM")
  maaslin_results_for_species$country <- factor(maaslin_results_for_species$country, levels = desired_order)
  maaslin_results_for_species$country <- as.character(maaslin_results_for_species$country)
  # add the 95% confidence intervals
  maaslin_results_for_species <- maaslin_results_for_species %>%
    mutate(
      ci_lower = coef - 1.96 * stderr,
      ci_upper = coef + 1.96 * stderr
    )
  return(maaslin_results_for_species)
}


get_prevalence <- function(strataa_metaphlan_data_species, groups_to_analyse) {
  patch_metaphlan_data <- read.csv(file = patch_metaphlan_handle, header= TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names=FALSE)
  patch_metaphlan_data$lowest_taxonomic_level <- sapply(str_split(row.names(patch_metaphlan_data), "\\|"), function(x) x[length(x)])
  patch_metaphlan_data_species <- patch_metaphlan_data %>% filter(str_starts(lowest_taxonomic_level, 's__'))

  patch_metaphlan_data_species_temp <- patch_metaphlan_data_species %>% select(!c(lowest_taxonomic_level))

  patch_metadata <- read.csv(file = patch_metadata_handle, header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
  patch_metadata$SampleID <- row.names(patch_metadata)
  patch_metadata <- patch_metadata %>% filter(Group %in% c('Typhi', 'Paratyphi'), day_group == 'baseline') %>% select(c(isolate_ID, Diagnosis)) %>% rename(SampleID = isolate_ID, Group = Diagnosis)
  patch_prevalence <- calculate_prevalence(patch_metaphlan_data_species_temp, patch_metadata, NA, c('Typhi', 'Paratyphi'), 'patch') %>% ungroup()
  # change disease in value column to Acute typhoid and no_disease to Household contact, so that they match the strataa data
  patch_prevalence$value <- ifelse(patch_prevalence$value == 'disease', 'Acute typhoid', patch_prevalence$value)
  patch_prevalence$value <- ifelse(patch_prevalence$value == 'no_disease', 'Household contact', patch_prevalence$value)
  patch_prevalence['country'] <- 'CHIM'

  strataa_metaphlan_data_species_temp <- strataa_metaphlan_data_species %>% select(!c(lowest_taxonomic_level))
  metadata_temp <- metadata %>% select(c(SampleID, Country, Group)) 
  row.names(metadata_temp) <- metadata_temp$SampleID
  # View(metadata_temp)
  # View(strataa_metaphlan_data_species_temp)
  bgd_prevalence <- calculate_prevalence(strataa_metaphlan_data_species_temp, metadata_temp, 'Bangladesh', groups_to_analyse, 'strataa') %>% ungroup()
  bgd_prevalence['country'] <- 'Bangladesh'
  variables_for_analysis <- c("Group", "Sex", "Age", "Antibiotics_taken_before_sampling_assumptions", "sequencing_lane")
  mwi_prevalence <- calculate_prevalence(strataa_metaphlan_data_species_temp, metadata_temp, 'Malawi', groups_to_analyse, 'strataa') %>% ungroup()
  mwi_prevalence['country'] <- 'Malawi'
  variables_for_analysis <- c("Group", "Sex", "Age", "Antibiotics_taken_before_sampling_assumptions")
  nep_prevalence <- calculate_prevalence(strataa_metaphlan_data_species_temp, metadata_temp, 'Nepal', groups_to_analyse, 'strataa') %>% ungroup()
  nep_prevalence['country'] <- 'Nepal'
  # View(bgd_prevalence)
  prevalence <- rbind(bgd_prevalence, mwi_prevalence, nep_prevalence, patch_prevalence)

  prevalence$lowest_taxonomic_level <- sapply(str_split(prevalence$feature, "\\|"), function(x) x[length(x)])
  return(prevalence)
}


forest_plot <- function(prevalence, species_of_interest, species_maaslin_results){
  # species_prevalence <- prevalence %>% filter(grepl(species_of_interest, lowest_taxonomic_level)) %>% select(!c(feature, variable))
  species_prevalence <- prevalence %>% filter(lowest_taxonomic_level %in% species_of_interest) %>% select(!c(feature, variable))
  # i want to pivot pcopri_species wider
  species_prevalence <- species_prevalence %>% pivot_wider(names_from = 'value', values_from = 'Median_Prevalence')

  species_maaslin_results <- species_maaslin_results %>% left_join(species_prevalence, by =c("lowest_taxonomic_level", 'country'))

  species_maaslin_results <- species_maaslin_results %>%
    mutate(across(c(coef, stderr, ci_lower, ci_upper, `Acute typhoid`, `Household contact`), ~ formatC(.x, format = "fg", digits = 3)))
  species_maaslin_results <- species_maaslin_results %>%
    mutate(across(c(pval, qval), ~ formatC(.x, format = "fg", digits = 2)))

  tabletext <- cbind(
    lowest_taxonomic_level = species_maaslin_results$lowest_taxonomic_level,
    country = species_maaslin_results$country,
    N = species_maaslin_results$N,
    coef = species_maaslin_results$coef,
    stderr = species_maaslin_results$stderr,
    qval = species_maaslin_results$qval,
    abundunce_disease = species_maaslin_results$`Acute typhoid`,
    abundunce_control = species_maaslin_results$`Household contact`
  )

  mean <- as.numeric(species_maaslin_results$coef)
  lower <- as.numeric(species_maaslin_results$coef) - 1.96 * as.numeric(species_maaslin_results$stderr)
  upper <- as.numeric(species_maaslin_results$coef) + 1.96 * as.numeric(species_maaslin_results$stderr)
  # remove NaN values from lower and upper, dont set them as 0, remove them entirely

  x_min <- min(na.omit(lower)) - 1  # Extend range slightly for better visualization
  x_max <- max(na.omit(upper)) + 1
  x_ticks <- seq(floor(x_min), ceiling(x_max), by = 2)  # Adjust 'by' for desired interval

  # Plot the forest plot with table
  forestplot(
    labeltext = tabletext,
    mean = mean,
    lower = lower,
    upper = upper,
    xlab = "Effect Size",
    xticks = x_ticks,
    new_page = TRUE,
    is.summary = c(FALSE, FALSE, FALSE, FALSE),
    txt_gp = fpTxtGp(label = gpar(cex = 0.8), ticks = gpar(cex = 0.7), xlab = gpar(cex = 0.8)),
    colgap = unit(0.2, "cm"),
    col = fpColors(box = "royalblue", lines = "darkblue", summary = "royalblue")) %>% 
      fp_set_zebra_style("#EFEFEF") %>% 
      fp_add_header(country = "Cohort" %>% fp_align_center(),
                  lowest_taxonomic_level = "Species" %>% fp_align_center(),
                  coef = 'Coefficient' %>% fp_align_center(),
                  stderr = 'Standard\nError' %>% fp_align_center(),
                  qval = 'Q-value' %>% fp_align_center(),
                  N = "N" %>% fp_align_center(),
                  abundunce_disease = "Median\nabundance\ntyphoid" %>% fp_align_center(),
                  abundunce_control = "Median\nabundance\ncontrol" %>% fp_align_center()
                  )
}


run_forest_plot <- function(prevalence, species_of_interest, combined_maaslins){
  species_maaslin_results <- get_maaslin_results_for_species(combined_maaslins, species_of_interest)
  forest_plot(prevalence, species_of_interest, species_maaslin_results)
}
