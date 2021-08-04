#this function reads the bracken output files and puts them all in a table
#Note: this code is dependent from the name of the output files (ie the naming convention of the samples). The regexp might needs modifying 
get_bracken_summarise <- function(folder, filename_regex){
  # compile all reports
  bracken_reports_files <- list.files(folder, full.names = T)
  
  combined_bracken <- NULL
  for (file in 1:length(bracken_reports_files)){
    file_info <- file.info(bracken_reports_files[file])
    # if the file contains data then add it to the combined table
    if (file_info$size > 99){
      bracken_reports_file <- read_tsv(bracken_reports_files[file])
      #browser()
      #here change the reg expr to match the sample names
      name <- str_extract(rownames(file_info), filename_regex)
      bracken_reports_file$sample_ID <- name
      combined_bracken <- rbind(combined_bracken, bracken_reports_file)
    }
  }
  return(combined_bracken)
}

#this function reads the summarised data, plots them, filters low abundunce reads, save and plot the filtered data
filter_data <- function(folder, prefix, filename_regex, output){
  #browser()
  
  #create the names of the parsed bracken files - one with all the data and one with the filtered
  summary_file <- paste(output, "summarised_", prefix, "_kraken.txt", sep = "")
  otu_file <- paste(output, "summarised_", prefix, "_otu.txt", sep = "")
  filtered_otu_file <- paste(output, "summarised_filtered_", prefix, "_otu.txt", sep = "")
  
  #if the filtered data exists, load it and exit
  if (file.exists(filtered_otu_file)){
    data <- read.csv(filtered_otu_file, header=T, sep = "\t")
    return(data)
  }
  
  #read the bracken output and put it in 1 file
  data <- get_bracken_summarise(folder, filename_regex)
  write.table(data, summary_file, row.names = F, quote=FALSE, sep='\t')
  
  
  #create the otu matrix
  otu <- data.frame(acast(data, name ~ sample_ID, value.var = "kraken_assigned_reads", fill=0))
  colnames(otu) = gsub(pattern = "X", replacement = "", x = colnames(otu))
  write.table(otu, otu_file, row.names = T, quote=FALSE, sep='\t')
  
  pdf(file = paste(output, prefix, "_abundances_hist.pdf", sep = ""))
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
    
  pdf(file = paste(output, prefix, "_filtered_abundances_hist.pdf", sep = ""))
  hist(as.matrix(otu_table_rare_removed_norm_cpm), 
       main="Histogram For filtered OTU table", 
       breaks = 100, 
       xlab="Relative Abundance", 
       border="blue")
  dev.off()
  
  
  return(filtered_otu)
}


bracken_pca <- function(data, meta, prefix, meta_column, out_folder, show_label = F){
  #browser()  
  #bracken_matrix_rev <- acast(data[1:8], sample_ID ~ name, value.var = "fraction_total_reads", fill=0)
  level <- str_split(out_folder, "_", simplify = TRUE)
  level <- str_remove(level[1,2],"/")
  
  out_folder <- paste(out_folder, "2_PCA/", sep = "")
  if (!dir.exists(out_folder)){ print("Creating folder: "); dir.create(out_folder) }
  
  filename <- paste(out_folder, prefix, " pca.pdf", sep = "")
  
  pca <- prcomp(as.data.frame(data), scale. = F)
  
  print(autoplot(pca, label = show_label, label.size = 3, scale = 0, data = meta, colour = meta_column) + ggtitle(paste("PCA","-", prefix,"(", level[1],  "level", ")")))
  ggsave(filename)
  
  return(pca)
  
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


clr_transform <- function(data, summary_column, filename, out_folder){

  #browser()
  #matrix <- acast(data, name ~ sample_ID, value.var = summary_column, fill=0)
  matrix <- data
  
  #this transposes the table to function properly, logs the values and replaces the 0s with a very low logarthmic number
  matrix.czm <- cmultRepl(t(matrix), label=0, method="CZM")
  
  # The table needs to be transposed again (samples as COLUMNS) - clr is calculated - output will have samples as rows
  matrix.clr <- t(apply(matrix.czm, 1, function(x){log(x) - mean(log(x))}))
  
  out_folder <- paste(out_folder, "1_normalized_data/", sep = "")
  if (!dir.exists(out_folder)){ dir.create(out_folder) }
  
  out_file <- paste(out_folder, filename, sep = "")
  write.table(matrix.clr, out_file, row.names=T, col.names=T, sep = "\t")
  
  return(matrix.clr)
}


tmm_transformation <- function(data, summary_column, filename, out_folder){
  #browser()
  #matrix <- acast(data, name ~ sample_ID, value.var = summary_column, fill=0)
  matrix <- data
  
  dgeObj <- DGEList(matrix)
  
  tmm <- calcNormFactors(dgeObj, method="TMM", na.action = na.omit)
  
  tmmScaleFactors <- dgeObj$samples$lib.size * tmm$samples$norm.factors
  tmmExp <- round(t(t(tmm$counts)/tmmScaleFactors) * mean(tmmScaleFactors))
  

  logcounts <- t(cpm(tmmExp,log=TRUE))
  
  out_folder <- paste(out_folder, "1_normalized_data/", sep = "")
  if (!dir.exists(out_folder)){ dir.create(out_folder) }
  
  out_file <- paste(out_folder, filename, sep = "")
  write.table(logcounts, out_file, row.names=T, col.names=T, sep = "\t")
  
  return(logcounts)
  
}


#this operates on a braken summary matrix
relative_abundances_histogram <- function(data, level, out_folder){
  
  est_reads_sum <- data %>% group_by(name) %>% summarise(total_reads = sum(new_est_reads)) %>% mutate(relative_abundance = total_reads / sum(total_reads) * 100)
  est_reads_sum <- transform(est_reads_sum, name = reorder(name, -relative_abundance))
  est_reads_sum <- est_reads_sum[order(-est_reads_sum$relative_abundance),]
  g <- ggplot(est_reads_sum[1:25,], aes(x=name, y=relative_abundance)) + geom_bar(stat="identity") + labs(x=level, y="Relative abundance") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g <- ggplotly(g)
  
  
  withr::with_dir(out_folder, saveWidget(g, file = "overall_abundances_histogram.html"))
  
  
  
}


calculate_beta <- function(data, meta, group, output_folder, prefix, level, print_extra_plots = FALSE){
  
  #browser()
  
  #make the first column headers
  rownames(meta) <- meta[,1]
  #meta <- meta %>% remove_rownames %>% column_to_rownames()
  
  #and order
  meta <- meta[ order(row.names(meta)), ]
  
  
  #order data table as well to be sure
  data <- data[ order(row.names(data)), ]
  rownames(data) <- gsub("#","_",rownames(data))

  #make sure that you have the same ids
  common.ids <- intersect(rownames(meta), rownames(data))
  
  #make sure you have the correct meta
  meta <- meta[common.ids,]
  data <- data[common.ids,]
  
  #transform the table to relative abundances
  data <- sweep(data, 1, rowSums(data),'/')
  
  #calculate bray curtis distance matrix
  d.bray <- vegdist(data)
  #transform it to a matrix to save it
  beta_matrix <- as.matrix(d.bray)
  
  #Perform PCoA
  pc.bray <- cmdscale(d.bray,k=2, eig = T)
  pcoa.var <- round(pc.bray$eig/sum(pc.bray$eig)*100, 1)
  pcoa.values <- pc.bray$points
  pcoa.data <- data.frame(Sample = rownames(pcoa.values), X=pcoa.values[,1], Y = pcoa.values[,2])
  #get the metadata column to paint the plot 
  paint <- meta[,group]
  #shapes <- meta[,9]
  
  
  
  output_folder <- paste("./", output_folder, "4_beta/", sep = "")
  if (!dir.exists(output_folder)){ dir.create(output_folder) }
  title <- paste("Beta diversity PCoA: ", level, " level", sep = "")
  
  #save the table
  out_file <- paste(output_folder, "pairwise_beta.txt", sep = "")
  write.table(beta_matrix, out_file, row.names=T, col.names=T, sep = "\t")
  
  
  #plot and save
  file_path <- paste(output_folder, prefix, "_beta_PCoA.pdf", sep = "")
  
  g1 <- ggplot(pcoa.data, aes(x=X, y=Y, colour = paint)) + 
          ggtitle(title) + xlab(paste("MDS1 - ", pcoa.var[1], "%", sep="")) + ylab(paste("MDS2 - ", pcoa.var[2], "%", sep="")) + 
          guides(colour=guide_legend(title=prefix)) +
          geom_point() #+ geom_text(aes(label=Sample),hjust=0, vjust=0)
  ggsave(file_path)
  
  print(g1)
  
}



calculate_alpha <- function(data, summary_column, meta, group, output_folder, prefix, level, print_extra_plots = FALSE){
  
  #browser()
  
  #calculate alpha diverities 
  #input_matrix <- acast(data, sample_ID ~ name, value.var = summary_column, fill=0)
  input_matrix <- t(data)
  
  alpha <- vegan::diversity(input_matrix, index="shannon")
  alpha_table <- tibble(isolate = names(alpha), alpha = alpha)
  alpha_table <- left_join(meta, alpha_table, by="isolate")
  
  
  #for t.test
  #all pairwise combinations
  my_comparisons <- combn(levels(meta[,eval(group)]), 2, simplify = F)
  
  #pairwise t test
  f <- paste("alpha~", group, sep = "")
  pv <- compare_means(as.formula(f),  data = alpha_table, method = "t.test")
  
  gr <- pv$p <= 0.05
  
  
  file_path <- paste(output_folder, "3_alpha/", prefix, "_alpha.pdf", sep = "")
  output_folder <- paste("./", output_folder, "3_alpha/", sep = "")
  if (!dir.exists(output_folder)){ dir.create(output_folder) }
  
  title <- paste("Alpha div.(T test): ", prefix, sep = "")
  
  if(all(gr == FALSE)){
    g1 <- ggplot(alpha_table, aes_(x=as.name(group), y=alpha)) + 
        ggtitle(title) + 
        labs(x="", y="Alpha diversity") +
        geom_boxplot(fill="slateblue", alpha=0.8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_point(alpha = 0.3, position = "jitter") 
     
  
    ggsave(file_path)
  } else {
    g1 <- ggplot(alpha_table, aes_(x=as.name(group), y=alpha)) + 
        ggtitle(title) + 
        labs(x="", y="Alpha diversity") +
        geom_boxplot(fill="slateblue", alpha=0.8) +
        geom_point(alpha = 0.3, position = "jitter") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_signif(comparisons = my_comparisons[gr], step_increase = 0.1, annotation = pv$p.signif[gr]) 
    
        ggsave(file_path)
  }
  
  #############################################
  
  #ks.test
  result_table <- combn(levels(meta[,eval(group)]), 2, simplify = F)
  
  result_table <- as.data.frame(result_table)
  result_table <- t(result_table)
  result_table <- as.data.frame(result_table)
  rownames(result_table) <- NULL
  result_table$D <- NA
  result_table$p_value <- NA
  result_table$significance <- NA
  
  
  for (row in 1:nrow(result_table)){
    
    x <- alpha_table[alpha_table[,group]== toString(result_table[row,1]),]
    y <- alpha_table[alpha_table[,group]== toString(result_table[row,2]),]
    
    
    result <- ks.test(x$alpha , y$alpha)
    result_table[row,3] <- result$statistic
    result_table[row,4] <- result$p.value
    
    
    if( result$p.value >= 0.05 ) {result_table[row,5] <- NA}
    if(result$p.value < 0.05 & result$p.value >= 0.01) {result_table[row,5] <- "*"}
    if(result$p.value < 0.01 & result$p.value >= 0.001) {result_table[row,5] <- "**"}
    if(result$p.value < 0.001 & result$p.value >= 0.0001) {result_table[row,5] <- "***"}
    if(result$p.value < 0.0001) {result_table[row,5] <- "****"}
  }
  
  
  my_comparisons <- combn(levels(meta[,eval(group)]), 2, simplify = F)
  gr <- result_table$p_value <= 0.05
  my_comparisons[gr]
  
  
  file_path <- paste(output_folder, prefix, "_alpha_ks.pdf", sep = "")
  
  
  title <- paste("Alpha div.(KS test): ", prefix, sep = "")
  
  
  if(all(gr == FALSE)){
  g2 <- ggplot(alpha_table, aes_(x=as.name(group), y=alpha)) + 
        ggtitle(title) + labs(x="", y="Alpha diversity") + 
        geom_boxplot(fill="slateblue", alpha=0.8) +
        geom_point(alpha = 0.3, position = "jitter") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file_path)
  }
  
  else{
    g2 <- ggplot(alpha_table, aes_(x=as.name(group), y=alpha)) + 
      ggtitle(title) + labs(x="", y="Alpha diversity") +
      geom_boxplot(fill="slateblue", alpha=0.8) +
      geom_point(alpha = 0.3, position = "jitter") +
      geom_signif(comparisons = my_comparisons[gr], step_increase = 0.1, annotation = result_table$significance[gr]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file_path)
    
  }
  #############################################  
  print(grid.arrange(g1, g2, ncol=2))
  
  
  
  #histogram
  pdf(file = paste(output_folder, "alpha_distribution.pdf", sep = ""))
  print(hist(alpha_table$alpha, main="Shannon diversity histogram", xlab="", breaks=30))
  dev.off()
  
  #accumulation curve
  richness <- specaccum(input_matrix, "random", permutations=100)
  pdf(file = paste(output_folder, "alpha_accumulation.pdf", sep = ""))
  plot(richness,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab="Number of Samples", ylab=paste(level, "Richness", sep=" "), main="Accumulation curve")
  #boxplot(richness, col="yellow", add=TRUE, pch="+")
  dev.off()
  
  if(print_extra_plots){
  print(hist(alpha_table$alpha, main="Shannon diversity histogram", xlab="", breaks=30))
    plot(richness,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlab="Number of Samples", ylab=paste(level, "Richness", sep=" "), main="Accumulation curve")
  }
  
 
}

create_bar_plot <- function(excel_input, column, out_folder){
  #browser()
 
  
  #load data
  otu_mat<- read_excel(excel_input, sheet = "unnormalised_OTU")
  tax_mat<- read_excel(excel_input, sheet = "TAX_map")
  samples_df <- read_excel(excel_input, sheet = "metadata")
  
  #reformat the tables a bit
  row.names_otu_mat <- otu_mat$OTU
  otu_mat <- otu_mat %>% dplyr::select(-1) 
  row.names(otu_mat) <- row.names_otu_mat
          
  row.names_tax_mat <- tax_mat$OTU
  tax_mat <- tax_mat %>% dplyr::select (-OTU) 
  row.names(tax_mat) <- row.names_tax_mat
  
  row.names_samples_df <- samples_df$ID
  samples_df <- samples_df %>% dplyr::select (-ID) 
  row.names(samples_df) <- row.names_samples_df 
  
  #turn them to the proper format
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  samples_df <- as.data.frame(samples_df)
  
  #add pseudocounts 
  otu_mat <- otu_mat[,] + 0.5
  
  #turn to pseq input format 
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df,  )
  
  #create pseq obj
  pseq <- phyloseq(OTU, TAX, samples)
  
  
  #Normalize number of reads in each sample using median sequencing depth.
  total = median(sample_sums(pseq))
  standf = function(x, t=total) round(t * (x / sum(x)))
  pseq = transform_sample_counts(pseq, standf)
  
  
  
  #bar plots  
  p1 <-plot_bar(pseq, fill="Class") + 
    facet_wrap(as.formula(paste("~", column)), scales="free_x", nrow=1) +
    geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack", show.legend = FALSE) +
    theme(axis.text.x=element_blank()) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 8)
    ) 
          
  print(p1)
  ggsave(paste(out_folder, column, "_barplot.pdf", sep = ""))
  
  #average and plot 
  pseq_fraction <- merge_samples(pseq, column)
  sample_data(pseq_fraction)$column <- levels(sample_data(pseq_fraction)$column)
  pseq_fraction = transform_sample_counts(pseq_fraction, function(x) 100 * x/sum(x))
  
  p2 <- plot_bar(pseq_fraction, fill = "Class") + 
    geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
  
  print(p2)
  ggsave(paste(out_folder, column, "_avg_barplot.pdf", sep = ""))
  
 
}
  
plot_glm <- function(glm_result, out_file, n, threshold, title, neg_log, pos_log, neg_color, pos_color){
  
  #browser()
  
  result <- glm_result
  
  #plot
  plotSmear(result, de.tags = rownames(result))
  
  resLRTfilt <- topTags(result, n = nrow(result$counts))
  volcanoData <- cbind(resLRTfilt$table$logFC, -log10(resLRTfilt$table$FDR), rownames(resLRTfilt$table))
  colnames(volcanoData) <- c("logFC", "negLogPval", "OTU")
  DEGs <- resLRTfilt$table$FDR < threshold & abs(resLRTfilt$table$logFC) > 1
  point.col <- ifelse(DEGs, "red", "black")
  plot(volcanoData, pch = 16, col = point.col, cex = 0.5)
  
  significant_only <- volcanoData[DEGs, ]
  significant_only <- head(significant_only,n)
  direction.col <- ifelse(significant_only[,1] > 0, pos_log, neg_log)
  significant_only <- data.frame(significant_only)
  significant_only$negLogPval <- as.numeric(as.matrix(significant_only$negLogPval))
  significant_only$Overexpressed <- direction.col
  
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  
  p<-significant_only %>% 
    mutate(OTU = fct_reorder(OTU, negLogPval)) %>% #use dplyr::desc to change
    ggplot(aes(x=negLogPval, y=OTU, col = Overexpressed  )) + geom_point(size=1) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    scale_x_continuous(labels=scaleFUN) + 
    scale_color_manual(values=c(neg_color, pos_color)) +
    #scale_color_manual(values=c( "magenta", "cyan")) +
    ylab("") + xlab("-log(Pvalue)") + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  plot(p)
  
  
  tree_input <- topTags(result, n=Inf)$table
  SIGs <- tree_input$FDR < threshold & abs(tree_input$logFC) > 1
  SIGs.colour <- ifelse(SIGs, "Significant", "Not Significant")
  tree_input <- cbind(gene=rownames(tree_input ), tree_input )
  
  p2 <- ggplot(data = tree_input, aes(x=logFC, y=-log10(PValue)))  +
          geom_point(aes(col=SIGs.colour)) + 
          scale_color_manual(values=c("black", "red")) +
          #ggtitle("Volcano Plot") + 
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_text_repel(data=head(tree_input, 15), aes(label=gene), size=2.7) + 
          theme_minimal()
  
  
  #plot in file
  pdf(file = out_file)
  
  #plotSmear(result, de.tags = rownames(result))
  #plot(volcanoData, pch = 16, col = point.col, cex = 0.5) 

  plot(p)
  plot(p2)
  
  dev.off()
  
}

#plot one OTU as a box plot
plot_one_gene <- function(glm_results, otu_table, meta_table, meta_column, filename){
  #browser()
  
  plist <- list()
  sig_taxa <- rownames(glm_results)
  i <- 1
  
  for (taxon in sig_taxa){
    print(taxon)
    
  
    one_otu <- as.data.frame(cpm(otu_table[, taxon], log = TRUE))
    one_otu$Diagnosis <- meta_table[rownames(one_otu), meta_column]
  
    plist[[i]] <- ggplot(one_otu, aes(x=Diagnosis, y=V1)) + 
      geom_boxplot(fill="slateblue", alpha=0.8) + 
      geom_point(alpha = 0.3, position = "jitter") +
      ggtitle(taxon) + theme(plot.title = element_text(hjust = 0.5)) +
      labs(y="Normalised Expression levels") +
      geom_signif(comparisons = list(c("0", "1")),map_signif_level=TRUE)
    
                
    i <- i+1;
    
  }
  
  #grid.arrange(grobs = plist, ncol = 2) ## display plot
  
  #pdf(file = "../1_species/5_glm/edgeR/baseline_week1pre/test.pdf")
  ggsave(filename, marrangeGrob(plist, nrow = 2, ncol = 2), width = 297, height =  210, units = "mm")
  
}