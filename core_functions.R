
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