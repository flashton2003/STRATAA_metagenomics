library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)



patch_input_data = read.csv(file = "~/Dropbox/GordonGroup/STRATAA_Microbiome/patch_data/metaphlan/results/2023.05.05/2023.05.05.patch.metaphlan.csv", header= TRUE, sep = ",", stringsAsFactors = FALSE, check.names=FALSE)


patch_metadata <- read.csv(file = "~/Dropbox/GordonGroup/STRATAA_Microbiome/patch_data/metaphlan/results/2023.05.05/samples_n383_allMetadata_SubsettedNewVariables_040318_KH.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#View(patch_input_data)
# View(patch_metadata)

# pivot longer so that every row has species, one sample and the rpkm
patch_input_data_lg <- patch_input_data %>% pivot_longer(!clade_name, names_to = "sample_name", values_to = "rpkm")

# do a left join to add the metadata to the patch_input_data_lg, only select the columns we want - time, group, diagnosis
patch_metadata_subset <- patch_metadata %>% select(isolate_ID, Patient, Day, day_group, Group, Diagnosis)
patch_input_data_lg <- left_join(patch_input_data_lg, patch_metadata_subset, by = c("sample_name" = "isolate_ID"))
# View(head(patch_input_data_lg, 10))

# filter so that the the clade name contains Prevotella_copri_clade_A
patch_prev_A <- patch_input_data_lg %>% filter(grepl("Prevotella_copri_clade_A$", clade_name)) %>% group_by(Patient) %>% filter(any(rpkm > 0))

patch_prev_A$Patient <- factor(patch_prev_A$Patient)

# filter patch_input_data_lg_prevotella_copri_clade_a so that at least one of the rpkm values is greater than 0.1 for each sample_name
#patch_has_prev_A <- patch_prev_A %>% group_by(Patient) %>% filter(any(rpkm > 0))

# filter patch_input_data_lg_prevotella_copri_clade_a so that we have sample_names that have at least two 

# test <- patch_input_data_lg_prevotella_copri_clade_a %>% filter(Patient == 8018)
# test$day_group <- factor(test$day_group, levels = c("baseline", "1-2 week", "one month", "three month"))



# filter patch_input_data_lg_prevotella_copri_clade_a to have only typhi patients
patch_prev_A <- arrange(patch_prev_A, Group)
patch_prev_A$Patient <- factor(patch_prev_A$Patient)

patch_prev_A_typhoid <- patch_prev_A %>% filter(Diagnosis == "typhi")
patch_prev_A_no_typhoid <- patch_prev_A %>% filter(Diagnosis == "no typhi")
patch_prev_A_bicarbonate <- patch_prev_A %>% filter(Group == 'Bicarbonate')





# draw a graph that has the rpkm on the y axis, the time on the x axis, and the color is the diagnosis for Prevotella_copri_clade_A
# set the alpha to 0.5

# use a consistent colour scheme between plots, colouring based on Group
myColors <- brewer.pal(length(unique(patch_prev_A$Group)),"Set1")
names(myColors) <- levels(factor(patch_prev_A$Group))
colScale <- scale_colour_manual(name = "Group", values = myColors)

plot_prev_rpkm <- function(patch_prev_A_typhi, patient_id){
    patient_prev_A_typhi <- patch_prev_A_typhi %>% filter(Patient == patient_id)
    ggplot(data = patient_prev_A_typhi, aes(x = Day, y = rpkm)) + 
        geom_point(aes(colour = Group), alpha = 0.5) +
        # guides(colour = FALSE) +
        ggtitle(paste("Patient", patient_id)) +
        colScale
        # + geom_line()# + labs(title = "Prevotella_copri_clade_A", x = "Day Group", y = "RPKM")
}
# View(myColors)
# run plot_prev_rpkm on every patient in patch_prev_A_typhi, capture the output and make into a patchwork
patch_prev_A_typhoid_patients <- patch_prev_A_typhoid %>% distinct(Patient, .keep_all = TRUE) %>% arrange(Group)
patchwork::wrap_plots(lapply(patch_prev_A_typhoid_patients$Patient, plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_typhoid)) + patchwork::plot_annotation('Developed typhoid') + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

patch_prev_A_no_typhi_patients <- patch_prev_A_no_typhi %>% distinct(Patient, .keep_all = TRUE) %>% arrange(Group)
patchwork::wrap_plots(lapply( patch_prev_A_no_typhi_patients$Patient, plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_no_typhi)) + patchwork::plot_annotation("Didn't develop typhoid") + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")


patchwork::wrap_plots(lapply(unique(patch_prev_A_bicarbonate$Patient), plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_bicarbonate)) + patchwork::plot_annotation("Bicarbonate") + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
