library(readr)
library(ggplot2)
library(patchwork)

all_sites_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_running_3/5_glm/all_sites.age_gender_country.results_all.edgeR.tsv'
bangladesh_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_dhaka/5_glm/dhaka.age_gender.results_all.edgeR.tsv'
malawi_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_blantyre/5_glm/blantyre.age_gender.results_all.edgeR.tsv'
nepal_dpt_handle <- '/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/phil_kathmandu/5_glm/kathmandu.age_gender.results_all.edgeR.tsv'

bangladesh_dpt <- read_delim(bangladesh_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
malawi_dpt <- read_delim(malawi_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
nepal_dpt <- read_delim(nepal_dpt_handle, delim = "\t", escape_double = FALSE,  trim_ws = TRUE)


do_qqplot <- function(edger_results, country){
  observed <- sort(edger_results$FDR)
  lobs <- -(log10(observed))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  df <- data.frame(lobs, lexp)
  p <- ggplot(df, aes(x = lexp, y = lobs)) + geom_point() + ggtitle(country) + geom_abline()
  return(p)
}

b <- do_qqplot(bangladesh_dpt, 'Bangladesh')
m <- do_qqplot(malawi_dpt, 'Malawi')
n <- do_qqplot(nepal_dpt, 'Nepal')

b / m / n

hist(bangladesh_dpt$FDR)
hist(malawi_dpt$FDR)
hist(nepal_dpt$FDR)
