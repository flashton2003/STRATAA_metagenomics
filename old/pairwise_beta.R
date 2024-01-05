pairwise_beta <- function (otu, meta, baseline, timepoints_column, group_column, out_folder){

  
  
#give metadata file
meta <- meta %>% remove_rownames %>% column_to_rownames(var="isolate")
meta$isolate_ID <- gsub('#', '_', meta$isolate_ID)
sampledata <- sample_data(meta)

#OTU
otu <- t(otu)
otu <- otu_table(otu, taxa_are_rows = TRUE)

pseq_obj <- phyloseq(otu, sampledata)


metadata <- meta(pseq_obj)
time_points <- levels(metadata[,timepoints_column]) #maybe try to order that


for(time_point in time_points){
  if(time_point == baseline){next}
  #browser()
  
  #get the subset of baseline vs every other time point
  pseq_subset <- prune_samples(meta(pseq_obj)[,timepoints_column] == baseline | meta(pseq_obj)[,timepoints_column] == time_point, pseq_obj)
  
  
  

  betas_list <- list()
  
  #get the groupings 
  phenotypic_groups <- as.character(unique(meta(pseq_subset)[[group_column]]))
  
  for (group in phenotypic_groups) {
    #browser()
    #group_members <- subset(meta(pseq_subset), Species == group)
    group_members <- meta(pseq_subset)[meta(pseq_subset)[[group_column]] == group, ]
      
    beta_div <- c()
    
    for (patient in group_members$Patient) {
      #browser()
      patient <- as.character(patient)
      # Pick the samples for this subject
      patient_samples <- subset(group_members, Patient == patient)
      # Check that the subject has two time points
      if (nrow(patient_samples) == 2) {
        print(patient)
        #browser()
        sample <- as.character(patient_samples$isolate)
        # Here with just two samples we can calculate the
        # beta diversity directly
        abundances_sample1 <- as.numeric(abundances(pseq_subset)[, sample[1]])
        abundances_sample2 <- as.numeric(abundances(pseq_subset)[, sample[[2]]])
        
        beta_div[[patient]] <- divergence(abundances_sample1,abundances_sample2,method = "bray")
      }
    }
    betas_list[[group]] <- beta_div
  }
  
  #browser()
  df <- as.data.frame(unlist(betas_list))
  s<- rownames(df)
  si<- as.data.frame(s)
  si<- separate(si, s, into = c('names','s'))
  df1<- bind_cols(df, si)
  rownames(df1)<- df1$s ; df1$s<- NULL
  
  colnames(df1)[1] <- "beta"
  
  my_comparisons <- combn(levels(as.factor(df1[,"names"])), 2, simplify = F)
  pv <- compare_means(beta ~ names,  data = df1, method = "wilcox.test")
  gr <- pv$p <= 0.05
  
  
  if(all(gr == FALSE)){
    p<- ggplot(df1, aes(x = names, y = beta)) + 
    geom_boxplot(fill="slateblue", alpha=0.7) + 
    geom_point(alpha = 0.6, position = "jitter") +
    ggtitle(time_point) + ylab('Beta Diversity') + xlab('')
    plot(p)
  }  else{
    p<- ggplot(df1, aes(x = names, y = beta)) + 
      geom_boxplot(fill="slateblue", alpha=0.7) + 
      geom_point(alpha = 0.6, position = "jitter") +
      ggtitle(time_point) + ylab('Beta Diversity') + xlab('') +
      geom_signif(comparisons = my_comparisons[gr], step_increase = 0.1, annotation = pv$p.signif[gr])
    plot(p)
  }

  ggsave(paste(out_folder,"4_beta/pairwise_beta/", time_point, ".pdf", sep = ""))
}



}