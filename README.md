# ASAP_human

## General information

This repository contains the scripts used to produce clinical and the microbiome analysis and figures presented in the manuscript : "Microbiome signature of Parkinson’s disease in healthy and genetically at-risk individuals".

## Files overview

### For clinical data

-   ASAP_Clinical_Microbiome_NatureMed.Rmd : script containing the analysis perform on clinical data

### For microbiome data

-   meteor_commands.txt : Meteo2 software commands to perform taxonomy and functional profiling both on gut and oral catalogs

-   merge_gut_oral.R : script containing the commands to merge gut and oral abundance tables

-   crocodeel_commands.txt : CroCoDeEL software commands to perform cross-contamination analysis on all samples

-   ASAPh_fecal_paper_fig_github.Rmd : script containing the analysis perform on microbiome data

-   script_function_ASAPh_paper.R : script containing the functions used for analysis of microbiome data

### Environment and packages

The version R used is 4.4.1

-   renv.lock : contains exact versions of R packages that are used for the study

-   .Rprofile : initializes the project-specific renv environment and ensures that R uses the isolated package library defined for this project.
