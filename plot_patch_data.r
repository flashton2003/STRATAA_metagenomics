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
patch_metadata_subset <- patch_metadata %>% select(isolate_ID, Patient, Day, day_group, Group, Diagnosis, pre_post_antibiotics, days_antibiotics)
patch_input_data_lg <- left_join(patch_input_data_lg, patch_metadata_subset, by = c("sample_name" = "isolate_ID"))
# View(head(patch_input_data_lg, 10))

patch_antibiotic_day <- patch_metadata_subset %>% mutate(antibiotic_delta = ifelse(days_antibiotics < 0, abs(days_antibiotics), NA)) %>% mutate(antibiotic_delta = ifelse(days_antibiotics > 0, days_antibiotics * -1, antibiotic_delta)) %>% mutate(antibiotic_delta = ifelse(days_antibiotics == 0, 0, antibiotic_delta))

patch_antibiotic_day <- patch_antibiotic_day %>% mutate(day_of_antibiotic = Day + antibiotic_delta) %>% distinct(Patient, day_of_antibiotic)

# this is an assert that there is only one day_of_antibiotic per patient - at the moment there is well done lorena!
stopifnot(length(patch_antibiotic_day$day_of_antibiotic) == length(unique(patch_antibiotic_day$Patient)))


# filter so that the the clade name contains Prevotella_copri_clade_A
patch_prev_A <- patch_input_data_lg %>% filter(grepl("Prevotella_copri_clade_A$", clade_name)) %>% group_by(Patient) %>% filter(any(rpkm > 0))



# filter patch_input_data_lg_prevotella_copri_clade_a to have only typhi patients
# patch_prev_A <- arrange(patch_prev_A, Group)
patch_prev_A$Patient <- factor(patch_prev_A$Patient)

patch_prev_A_typhoid <- patch_prev_A %>% filter(Diagnosis == "typhi")
patch_prev_A_no_typhoid <- patch_prev_A %>% filter(Diagnosis == "no typhi")
patch_prev_A_bicarbonate <- patch_prev_A %>% filter(Group == 'Bicarbonate')


# going to figure out which day each person had their antibiotics on
# take the day of each row, and if the days_antibiotics is negative, then add it to the day, if the days_antibiotics is positive, subtract it from the days




# draw a graph that has the rpkm on the y axis, the time on the x axis, and the color is the diagnosis for Prevotella_copri_clade_A
# set the alpha to 0.5

# use a consistent colour scheme between plots, colouring based on Group
myColors <- brewer.pal(length(unique(patch_prev_A$Group)),"Set1")
names(myColors) <- levels(factor(patch_prev_A$Group))
colScale <- scale_colour_manual(name = "Group", values = myColors)

plot_prev_rpkm <- function(patch_prev_A_typhi, patient_id, patch_antibiotic_day){
    patient_prev_A_typhi <- patch_prev_A_typhi %>% filter(Patient == patient_id)
    antibiotic_day <- patch_antibiotic_day %>% filter(Patient == patient_id) %>% pull(day_of_antibiotic)
    ggplot(data = patient_prev_A_typhi, aes(x = Day, y = rpkm)) + 
        geom_point(aes(colour = Group, shape = pre_post_antibiotics)) +
        # guides(colour = FALSE) +
        ggtitle(paste("Patient", patient_id)) +
        colScale +
        theme(axis.title = element_blank()) +
        geom_vline(xintercept = antibiotic_day, linetype="dotted", 
                color = "blue")
        # scale_shape_manual(name = "pre_post_antibiotics", values = myShapes)
        # + geom_line()# + labs(title = "Prevotella_copri_clade_A", x = "Day Group", y = "RPKM")
}
# View(myColors)
# run plot_prev_rpkm on every patient in patch_prev_A_typhi, capture the output and make into a patchwork
patch_prev_A_typhoid_patients <- patch_prev_A_typhoid %>% distinct(Patient, .keep_all = TRUE) %>% arrange(Group)
patchwork::wrap_plots(lapply(patch_prev_A_typhoid_patients$Patient, plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_typhoid, patch_antibiotic_day = patch_antibiotic_day)) + patchwork::plot_annotation('Developed typhoid') + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

patch_prev_A_no_typhoid_patients <- patch_prev_A_no_typhoid %>% distinct(Patient, .keep_all = TRUE) %>% arrange(Group)
patchwork::wrap_plots(lapply( patch_prev_A_no_typhoid_patients$Patient, plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_no_typhoid, patch_antibiotic_day = patch_antibiotic_day)) + 
    patchwork::plot_annotation("Didn't develop typhoid") + 
    patchwork::plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")


patchwork::wrap_plots(lapply(unique(patch_prev_A_bicarbonate$Patient), plot_prev_rpkm, patch_prev_A_typhi = patch_prev_A_bicarbonate, patch_antibiotic_day = patch_antibiotic_day)) + patchwork::plot_annotation("Bicarbonate") + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

