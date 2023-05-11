library(dplyr)
library(edgeR)
library(ggplot2)
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


#source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/config.R")

read_metadata <- function(path_to_metadata){
  meta <- read.csv(path_to_metadata, header=T, sep = "\t")
  # i dont think this line does anything, probably just there for historic reasons
  names(meta)[names(meta) == "sample_ID"] <- "isolate"
  names(meta)[names(meta) == "ID"] <- "isolate"
  meta <- meta %>% mutate(age_bracket=cut(Age, breaks=c(0, 1, 5, 15, Inf), labels=c("0-1", "1-5", "6-15", ">15")))
  meta$group_country <- paste(meta$Group, meta$Country, sep = '_')
  meta$group_antibiotic <- paste(meta$Group, meta$Antibiotics_taken_before_sampling_yes_no_assumptions, sep = '_')
  
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
    s <- summary(aov(alpha ~ Country * Sex *Group * age_bracket * Antibiotics_taken_before_sampling_yes_no_assumptions, data = meta))
  } else {
    s <- summary(aov(alpha ~ Sex *Group * age_bracket * Antibiotics_taken_before_sampling_yes_no_assumptions, data = meta))
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




combine_and_compare_edgeRs <- function(to_combine, covar_initials, output_folder_all_three, output_folder_blantyre, output_folder_dhaka, output_folder_kathmandu, combined_output_root, the_date, comp) {
  # this function is atrociously repetitive.
  # should just merge all of them to start with and then do filters on multiple columns.
  # dpt is differentially present taxa
  handles <- sapply(to_combine, make_name, covars = covars, comp = comp)

  all_sites_dpt <- read_delim(handles[[output_folder_all_three]], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  all_sites_dpt <- rename(all_sites_dpt, c(species=...1, all_sites_logFC = logFC, all_sites_logCPM = logCPM, all_sites_LR = LR, all_sites_PValue = PValue, all_sites_FDR = FDR))
  
  bangladesh_dpt <- read_delim(handles[[output_folder_dhaka]], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  bangladesh_dpt <- rename(bangladesh_dpt, c(species=...1, bangladesh_logFC = logFC, bangladesh_logCPM = logCPM, bangladesh_LR = LR, bangladesh_PValue = PValue, bangladesh_FDR = FDR))
  bangladesh_dpt_sig <- filter(bangladesh_dpt, bangladesh_FDR <= 0.01)
  print('Bangladesh FDR <= 0.01')
  print(length(bangladesh_dpt_sig$species))
  bangladesh_dpt_sig_up <- filter(bangladesh_dpt_sig, bangladesh_logFC >= 1)
  bangladesh_dpt_sig_down <- filter(bangladesh_dpt_sig, bangladesh_logFC <= -1)
  bangladesh_dpt_sig_up$bangladesh_FDR <- formatC(bangladesh_dpt_sig_up$bangladesh_FDR, format = 'e', digits = 1)
  bangladesh_dpt_sig_down$bangladesh_FDR <- formatC(bangladesh_dpt_sig_down$bangladesh_FDR, format = 'e', digits = 1)
  
  malawi_dpt <- read_delim(handles[[output_folder_blantyre]], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  malawi_dpt <- rename(malawi_dpt, c(species=...1, malawi_logFC = logFC, malawi_logCPM = logCPM, malawi_LR = LR, malawi_PValue = PValue, malawi_FDR = FDR))  
  malawi_dpt_sig <- filter(malawi_dpt, malawi_FDR <= 0.01)
  print('Malawi FDR <= 0.01')
  print(length(malawi_dpt_sig$species))
  malawi_dpt_sig_up <- filter(malawi_dpt_sig, malawi_logFC >= 1)
  malawi_dpt_sig_down <- filter(malawi_dpt_sig, malawi_logFC <= -1)
  malawi_dpt_sig_up$malawi_FDR <- formatC(malawi_dpt_sig_up$malawi_FDR, format = 'e', digits = 1)
  malawi_dpt_sig_down$malawi_FDR <- formatC(malawi_dpt_sig_down$malawi_FDR, format = 'e', digits = 1)
  
  nepal_dpt <- read_delim(handles[[output_folder_kathmandu]], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  nepal_dpt <- rename(nepal_dpt, c(species=...1, nepal_logFC = logFC, nepal_logCPM = logCPM, nepal_LR = LR, nepal_PValue = PValue, nepal_FDR = FDR))  
  nepal_dpt_sig <- filter(nepal_dpt, nepal_FDR <= 0.01)
  print('Nepal FDR <= 0.01')
  print(length(nepal_dpt_sig$species))
  nepal_dpt_sig_up <- filter(nepal_dpt_sig, nepal_logFC >= 1)
  nepal_dpt_sig_down <- filter(nepal_dpt_sig, nepal_logFC <= -1)
  
  malawi_bangladesh_dpt <- read_delim(handles[[output_folder_blantyre_dhaka]], delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
  malawi_bangladesh_dpt <- rename(malawi_bangladesh_dpt, c(species=...1, malawi_bangladesh_logFC = logFC, malawi_bangladesh_logCPM = logCPM, malawi_bangladesh_LR = LR, malawi_bangladesh_PValue = PValue, malawi_bangladesh_FDR = FDR))  
  malawi_bangladesh_dpt_sig <- filter(malawi_bangladesh_dpt, malawi_bangladesh_FDR <= 0.01)
  print('malawi_bangladesh FDR <= 0.01')
  print(length(malawi_bangladesh_dpt_sig$species))
  malawi_bangladesh_dpt_sig_up <- filter(malawi_bangladesh_dpt_sig, malawi_bangladesh_logFC >= 1)
  malawi_bangladesh_dpt_sig_down <- filter(malawi_bangladesh_dpt_sig, malawi_bangladesh_logFC <= -1)
  malawi_bangladesh_dpt_sig_up$malawi_bangladesh_FDR <- formatC(malawi_bangladesh_dpt_sig_up$malawi_bangladesh_FDR, format = 'e', digits = 1)
  malawi_bangladesh_dpt_sig_down$malawi_bangladesh_FDR <- formatC(malawi_bangladesh_dpt_sig_down$malawi_bangladesh_FDR, format = 'e', digits = 1)
  
  combined_dpt <- left_join(all_sites_dpt, bangladesh_dpt, by = "species")
  combined_dpt <- left_join(combined_dpt, malawi_dpt, by = "species")
  combined_dpt <- left_join(combined_dpt, nepal_dpt, by = "species")
  combined_dpt <- left_join(combined_dpt, malawi_bangladesh_dpt, by = "species")
  
  combined_output_folder <- file.path(combined_output_root, 'dge')
  if (!dir.exists(combined_output_folder)){ dir.create(combined_output_folder) }
  
  write_csv(combined_dpt, file= file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'combined_dges.csv', sep = '.')))
  # do some inner joins of malawi and bangladesh to get the taxa that are dp between them
  # then do a left join to pull in the info on those taxa from teh combined malawi_bangladesh analysis
  # i.e. the edgeR that pooled samples from both sites.
  View(malawi_bangladesh_dpt)
  malawi_and_bangladesh_sig_up <- inner_join(bangladesh_dpt_sig_up, malawi_dpt_sig_up, by = 'species')
  malawi_and_bangladesh_sig_up <- malawi_and_bangladesh_sig_up %>% left_join(malawi_bangladesh_dpt, by = 'species')
  malawi_and_bangladesh_sig_down <- inner_join(bangladesh_dpt_sig_down, malawi_dpt_sig_down, by = 'species')
  malawi_and_bangladesh_sig_down <- malawi_and_bangladesh_sig_down %>% left_join(malawi_bangladesh_dpt, by = 'species')
  
  View(malawi_and_bangladesh_sig_up)
  
  
  write_csv(malawi_and_bangladesh_sig_up, file= file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'bang_mal_up.filtered.csv', sep = '.')))
  write_csv(malawi_and_bangladesh_sig_down, file= file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'bang_mal_down.filtered.csv', sep = '.')))
  
  print(paste('Jaccard of unfiltered taxa of Bang, Mal, Nep = ', jaccard(bangladesh_dpt$species, malawi_dpt$species, nepal_dpt$species)), sep = '')
  
  #   venn.diagram(x = list(bangladesh_dpt_sig_up$species, malawi_dpt_sig_up$species), category.names = c('bangladesh', 'malawi'), filename = '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/combined_dge/results/2022.06.07/2022.06.07.venn_diagram.fdr_0.01.upreg.png', euler.d = FALSE, scaled = FALSE)
  # unfiltered
  venn.diagram(x = list(bangladesh_dpt$species, malawi_dpt$species, nepal_dpt$species), category.names = c('Bangladesh', 'Malawi', 'Nepal'), filename = file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'venn_diagram.no_filter.png', sep = '.')), euler.d = FALSE, scaled = FALSE, height=2200, width=2200)
  # fdr 0.01
  venn.diagram(x = list(bangladesh_dpt_sig$species, malawi_dpt_sig$species, nepal_dpt_sig$species), category.names = c('Bangladesh', 'Malawi', 'Nepal'), filename = file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'venn_diagram.fdr_0.01.png', sep = '.')), euler.d = FALSE, scaled = FALSE, height=2200, width=2200)
  # fdr 0.01 and logfc > 1
  venn.diagram(x = list(bangladesh_dpt_sig_up$species, malawi_dpt_sig_up$species, nepal_dpt_sig_up$species), category.names = c('Bangladesh', 'Malawi', 'Nepal'), filename = file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'venn_diagram.fdr_0.01.upreg.png', sep = '.')), euler.d = FALSE, scaled = FALSE)
  # fdr 0.01 and logfc < -1
  venn.diagram(x = list(bangladesh_dpt_sig_down$species, malawi_dpt_sig_down$species, nepal_dpt_sig_down$species), category.names = c('Bangladesh', 'Malawi', 'Nepal'), filename = file.path(combined_output_folder, paste(the_date, comp, covar_initials, 'venn_diagram.fdr_0.01.downreg.png', sep = '.')), euler.d = FALSE, scaled = FALSE)
  
  output <- list(sig_up_cases_controls = malawi_and_bangladesh_sig_up, sig_down_cases_controls = malawi_and_bangladesh_sig_down)
  return(output)
}


new_combine_and_compare_edgeRs <- function(){
  
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


run_maaslin <- function(feature_data, metadata, output_root, country, groups_for_analysis, variables_for_analysis, norm, trans, reference_groups){
  ifelse(!dir.exists(output_root), dir.create(output_root), FALSE)
  metadata_to_analyse <- metadata %>% filter(Country == country, Group %in% groups_for_analysis)
  View(metadata_to_analyse)
  vars_for_dirname <- paste(variables_for_analysis, collapse = '.')
  output_dir <- file.path(output_root, paste(country, paste(groups_for_analysis, collapse = '_vs_'), vars_for_dirname, sep = '_'))
  Maaslin2(input_data = feature_data, input_metadata = metadata_to_analyse, analysis_method = "LM", min_prevalence = 0,
           normalization  = norm,
           transform = trans,
           output         = output_dir, 
           fixed_effects  = variables_for_analysis,
           reference = reference_groups)
}





combine_maaslins <- function(bangladesh_maaslin, malawi_maaslin){
  # thanks chatgpt!
  combined_df <- bangladesh_maaslin %>%
    inner_join(malawi_maaslin, by = c("feature", 'metadata', 'value'), suffix = c("_bang", "_mal"))
  
  # Filter the combined data frame based on the conditions for coef > 0
  filtered_df_positive_coef <- combined_df %>%
    filter(
      qval_bang < 0.05 & qval_mal < 0.05 & 
        (coef_bang > 0 & coef_mal > 0)
    )
  
  # Filter the combined data frame based on the conditions for coef < 0
  filtered_df_negative_coef <- combined_df %>%
    filter(
      qval_bang < 0.05 & qval_mal < 0.05 & 
        (coef_bang < 0 & coef_mal < 0)
    )
  
  # Return the two filtered data frames
  return(list(positive_coef = filtered_df_positive_coef, negative_coef = filtered_df_negative_coef, all_features = combined_df))
}















