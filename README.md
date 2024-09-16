# STRATAA metagenome analysis scripts

Authors - Leonardos Mageiros & Philip Ashton

The code here was used for the analysis of microbiomes for Ashton et al., 2024 manuscript availabile [here]([url](https://www.medrxiv.org/content/10.1101/2024.09.02.24312347v1)).

The code is not intended to work smoothly on other machines, as it contains some hardcoded paths and the documentation is somewhat lacking.

The html files are a good place to start if you're interested in exploring this code, as the graphs etc are integrated, and the code has been run with all dependencies and raw datafiles. You will need to download the html file and open it with a web browser.

# Code structure

00.core_functions - encodes many functions required for analysis, these are imported into the other scripts for running.

01.metadata_analysis - analysis of basic information like age, sex breakdown.

02.run_maaslin - this script imports the metaphlan and big-map results (which are run separately) and runs maaslin on them.

03.taxonomic_analysis - makes graphs, tables, etc for taxonomic analysis.

04.functional_analysis - makes graphs, tables, etc for functional analysis.

05.random_forest - does the random forest analysis.

06.amr_analysis - analysis of AMR genes.

07.patch_analysis - analysis of controlled human infection model samples.
