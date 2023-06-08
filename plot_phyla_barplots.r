library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)

source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/core_functions.R")
source("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/from_Leo/Leonardos_analysis/bin/config.R")


data <- read_csv("/Users/flashton/Dropbox/GordonGroup/STRATAA_Microbiome/metaphlan/results/2023.05.11/2023.05.11.all_strataa_metaphlan.csv")

phyla <- data %>% mutate(phylum = str_extract(clade_name, "p__[A-Za-z0-9]*"))
metadata <- read_metadata(metadata_handle)

phyla <- data %>% filter(grepl("p__", clade_name)) %>% filter(!grepl("c__", clade_name)) %>% pivot_longer(!clade_name, names_to = "sample", values_to = "relative_abundance")

# show that phyla counts sum to 100
# phyla_summed_by_sample <- phyla %>% group_by(sample) %>% summarise(relative_abundance = sum(relative_abundance))

phyla_to_exclude <- phyla %>% group_by(clade_name) %>% 
    summarise(count = sum(relative_abundance > 1)) %>% 
    filter(count < 31) %>% 
    pull(clade_name)

phyla_clean <- phyla %>%
  filter(!(clade_name %in% phyla_to_exclude))

metadata_select <- metadata %>% select(Lane, Group, Country)

phyla_clean_metadata <- phyla_clean %>% left_join(metadata_select, by = c("sample" = "Lane"))

group_order <- c("Acute_Typhi", "Control_HealthySerosurvey", "Carrier")




phyla_clean_nepal <- phyla_clean_metadata %>% filter(Country == "Nepal")

# thanks chatgpt!
phyla_clean_nepal_fct <- phyla_clean_nepal %>%
  mutate(Group = factor(Group, levels = group_order),  # Convert group to a factor with the desired order
         group_order_numeric = as.numeric(Group),  # Create a new numeric variable based on the order of group
         sample = fct_reorder(sample, group_order_numeric))  # Reorder sample based on group_order_numeric




phyla_clean_fct <- phyla_clean_metadata %>%
  mutate(Group = factor(Group, levels = group_order),  # Convert group to a factor with the desired order
         group_order_numeric = as.numeric(Group),  # Create a new numeric variable based on the order of group
         sample = fct_reorder(sample, group_order_numeric))  # Reorder sample based on group_order_numeric



plot_per_country_abundance <- function(phyla_clean_metadata, country){
    # thanks chatgpt!
    phyla_clean_country <- phyla_clean_metadata %>% filter(Country == country)

    phyla_clean_country_fct <- phyla_clean_country %>%
    mutate(Group = factor(Group, levels = group_order),  # Convert group to a factor with the desired order
            group_order_numeric = as.numeric(Group),  # Create a new numeric variable based on the order of group
            sample = fct_reorder(sample, group_order_numeric))  # Reorder sample based on group_order_numeric

    p <- ggplot(data = phyla_clean_country_fct, aes(x = sample, y = relative_abundance, fill = clade_name)) + 
        geom_bar(stat = "identity") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_x_discrete(breaks=phyla_clean_country_fct$sample, labels=phyla_clean_country_fct$Group) + 
        ggtitle(country) +
        ylim(0, 100)
    p
}


plot_per_country_abundance(phyla_clean_metadata, "Nepal")
plot_per_country_abundance(phyla_clean_metadata, "Bangladesh")
plot_per_country_abundance(phyla_clean_metadata, "Malawi")
