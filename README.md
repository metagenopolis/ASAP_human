# ASAP_human

## General information

This repository contains the scripts used to produce clinical and the microbiome analysis and figures presented in the manuscript : "Microbiome signature of Parkinson’s disease in healthy and genetically at-risk individuals".

## Files overview

### For clinical data

-   ASAP_SCHAPIRA_Clinical_Microbiome_NatureMed_2026_Clinical_code.Rmd : script containing the analysis perform on clinical data

-   ASAP_SCHAPIRA_Codebook_Clinical_Dataset : file containing the definition of the clinical variables used in the study

### For microbiome data

-   meteor_commands.txt : Meteor2 software commands to perform taxonomy and functional profiling both on gut and oral catalogs

-   crocodeel_commands.txt : CroCoDeEL software commands to perform cross-contamination analysis on all samples

-   prepare_ASAPh_fecal_data.R : script containing the commands to prepare data for microbiome analysis, including the quality control of samples, the computation of alpha-diversity indices, the merge of gut and oral abundance tables, the preparation of functional data

-   MSP_set_ref_oss_gss_status_612_ref_paper_20220622.RData : contains the taxonomic reference of species

-   ASAPh_fecal_paper_fig_github.Rmd : script containing the analysis performed on microbiome data

-   script_function_ASAPh_paper.R : script containing the functions used for analysis of microbiome data

### Environment, packages and datasets

The version R used is 4.4.1

-   KRT_Menozzi_ASAP_Clinical_Microbiome_March2026_complete.csv : contains the description of the datasets/protocols/softwares used in the study and their sources

-   renv.lock : contains exact versions of R packages that are used for the study

-   .Rprofile : initializes the project-specific renv environment and ensures that R uses the isolated package library defined for this project.
