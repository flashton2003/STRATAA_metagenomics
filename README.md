# STRATAA metagenome analysis scripts

Authors - Leonardos Mageiros (primary author) & Philip Ashton

## getting started

1. Clone this repo i.e. `git clone git@github.com:flashton2003/STRATAA_metagenomics.git`
2. Install all the required packages
	a. They can be found in the `STRATAA_analysis_species.Rmd` file.
3. Change the location of the files which are sourced in lines 61-63 to the location to which you downloaded them.
4. Change L76 to set the root.dir as the directory in which you saved the `STRATAA_analysis_species.Rmd` file.
5. Change L80 to define `meta` based on the location of the `full_meta.txt` file from google drive.
6. Change L91 to set the `bracken_folder` variable to be the location of the bracken output from the google drive.
7. Change L93  to set the `out_folder` variable to be the directory where you want to save the output.
