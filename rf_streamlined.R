library(broom)
library(tidyverse)
library(stringr)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(heatmaply)
library(ggplot2)
library(readr)
library(ggfortify)
library(cluster)
library(autoplotly)
library(htmlwidgets)
library(compositions)
library(zCompositions)
library(vegan)
library(dendextend)
library(dplyr)
library(ggpubr)
library(edgeR)
library(rmarkdown)
library(gridExtra)
library(microbiome)
library(phyloseq)
library(apeglm)
library(vsn)
library(DESeq2)
library(edgeR)
library(pROC)
#libraries for RF
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 

#this is a file with some functions we are going to use to anlyse our results
source("parse_bracken_functions.r")

meta <- read.csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", header=T, sep = "\t")
names(meta)[names(meta) == "sample_ID"] <- "isolate"
names(meta)[names(meta) == "ID"] <- "isolate"



run_rf <- function(braken_folder, level, out_folder, place) {
  #give filename patern as a reg exp
  samples_regexp = "\\d+_\\d+_\\d+" #for patch and strata
  ifelse(!dir.exists(out_folder), dir.create(out_folder), FALSE)
  ifelse(!dir.exists(file.path(out_folder, '1_species')), dir.create(file.path(out_folder, '1_species')), FALSE)
  
  data_table <- filter_data(braken_folder, level, samples_regexp, out_folder)
  
  #normalise the PATCH data
  #data_table.normalised.tmm <- tmm_transformation(data_table, "kraken_assigned_reads", "tmm_normalised_data.txt", out_folder)
  #rownames(data_table.normalised.tmm) = gsub(pattern = "X", replacement = "", x = rownames(data_table.normalised.tmm))
  data_table.normalised.clr <- t(clr_transform(data_table, "kraken_assigned_reads", "clr_normalised_data.txt", out_folder))
  colnames(data_table.normalised.clr) = gsub(pattern = "X", replacement = "", x = colnames(data_table.normalised.clr)) 
  
  ## RF
  #https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial
  
  metadata <- meta %>% remove_rownames %>% column_to_rownames(var="isolate")
  metadata <- metadata[order(row.names(metadata)), ]
  #add at the end the diagnosis
  otu_table_scaled_state <- data.frame(data_table.normalised.clr)
  colnames(otu_table_scaled_state) = gsub(pattern = "X", replacement = "", x = colnames(otu_table_scaled_state))
  otu_table_scaled_state <- as.data.frame(t(otu_table_scaled_state))
  
  otu_table_scaled_state$Diagnosis <- metadata[rownames(otu_table_scaled_state), "Disease"]
  otu_table_scaled_state$Group <- metadata[rownames(otu_table_scaled_state), "Group"]
  
  set.seed(151)
  
  x <- otu_table_scaled_state[otu_table_scaled_state$Diagnosis != "na", ]
  
  otu_table_scaled_state_no_placebo <- droplevels(otu_table_scaled_state[otu_table_scaled_state$Diagnosis != "na", ])
  #otu_table_scaled_state_no_placebo$Diagnosis <- as.numeric(otu_table_scaled_state_no_placebo$Diagnosis)
  otu_table_scaled_state_no_placebo$Diagnosis <- as.factor(otu_table_scaled_state_no_placebo$Diagnosis)
  
  otu_table_scaled_state_no_placebo$Group <- NULL
  otu_table_scaled_state_no_placebo <- otu_table_scaled_state_no_placebo[grepl("NA", rownames(otu_table_scaled_state_no_placebo))==F,]
  
  
  # run the model
  RF_state_classify <- randomForest( x=otu_table_scaled_state_no_placebo[,1:(ncol(otu_table_scaled_state_no_placebo)-1)] , y=otu_table_scaled_state_no_placebo[ , ncol(otu_table_scaled_state_no_placebo)] , ntree=5001, importance=TRUE, proximities=TRUE )
  print(RF_state_classify)
  RF_state_classify_imp <- as.data.frame( RF_state_classify$importance )
  RF_state_classify_imp$features <- rownames( RF_state_classify_imp )
  RF_state_classify_imp_sorted <- arrange( RF_state_classify_imp  , desc(MeanDecreaseAccuracy)  )
  RF_state_classify_imp_sorted_percentile <- mutate(RF_state_classify_imp_sorted, percentile = ntile(desc(RF_state_classify_imp_sorted$MeanDecreaseAccuracy), 100))
  
  png(filename = paste(out_folder, '2_RF/', place, '_MeanDecreaseAccuracy.png', sep = ''))
  barplot(RF_state_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")
  dev.off()
  
  
  
  
  ifelse(!dir.exists(file.path(out_folder, '2_RF')), dir.create(file.path(out_folder, '2_RF')), FALSE)
  
  write_csv(RF_state_classify_imp_sorted_percentile, paste(out_folder, '2_RF/', place, '_features.importance.csv', sep = ''))
  # print the curve for the model
  metadata2 <- metadata
  metadata2$sample_id <- rownames(metadata2)
  metadata2 <- metadata2[(metadata2$sample_id %in% rownames(RF_state_classify$votes)), ]
  
  rf.roc<-roc(metadata2$Disease[metadata2$Disease != "na"] , RF_state_classify$votes[,2])
  label <- auc(rf.roc)
  ggroc(rf.roc, lwd=0.8, col="blue")+
    geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed", lwd=0.5) +
    ggtitle("ROC curve - All samples") + theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x=0.15, y=0.2, label= label) 
  
  ggsave(paste(out_folder, "2_RF/", place, "_all_roc.pdf", sep = "")) 
  
  RF_state_classify_sig <- rf.significance( x=RF_state_classify ,  xdata=otu_table_scaled_state[,1:(ncol(otu_table_scaled_state)-1)] , nperm=100 , ntree=501 )
  #print(RF_state_classify_sig)
  rf_output <- list(rf = RF_state_classify, rf_sig = RF_state_classify_sig)
  #rf_output <- list(rf = RF_state_classify)
  return(rf_output)
}




#define the braken output folder, the level we are operating and a folder for all the output
all_braken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output/species"
all_out_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/random_forest/"
all_rf_state_classify <- run_rf(all_braken_folder, level, all_out_folder, 'all')

level = "species"
blantyre_dhaka_bracken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_blantyre_dhaka/species"
blantyre_dhaka_out_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/random_forest_blantyre_dhaka/"
blantyre_dhaka_rf_output <- run_rf(blantyre_dhaka_bracken_folder, level, blantyre_dhaka_out_folder, 'blantyre_dhaka')

blantyre_dhaka_rf_output$rf_sig

blantyre_bracken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_blantyre/species/"
blantyre_out_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/random_forest_blantyre/"
blantyre_output <- run_rf(blantyre_bracken_folder, level, blantyre_out_folder, 'blantyre')

dhaka_bracken_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/1_taxonomic_profiling/bracken_output_dhaka/species/"
dhaka_out_folder = "/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/random_forest_dhaka/"
dhaka_output <- run_rf(dhaka_bracken_folder, level, dhaka_out_folder, 'dhaka')

