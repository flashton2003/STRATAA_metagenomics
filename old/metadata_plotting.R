library(readr)
library(dplyr)
library(ggplot2)


full_meta <- read_delim("~/Dropbox/GordonGroup/STRATAA_Microbiome/from_leo/Leonardos_analysis/0_metadata/full_meta.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)


ggplot(full_meta, aes(x = Country, y = Age, fill = Group)) + geom_boxplot()
