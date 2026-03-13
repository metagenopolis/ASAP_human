############################################################################################
### Data management and preparation for the fecal baseline analysis with 464 individuals ###
############################################################################################
# Authors : Victoria Meslier, Yani Ren

library(dplyr)
library(vegan)
library(data.table)

###############################
# Load metadata  
###############################
# v12 - October 2025 - contains the clinical data for ASAPh individuals, build by Elisa Menozzi 
metadata_last <- xlsx::read.xlsx("./data/ASAP_human_metadata_file_559ind_20250915_v12.xlsx", sheetIndex = 1, stringsAsFactors = FALSE) 
dim(metadata_last) # 559 individuals ; 129 variables 

# replace "NA" into real NA in each column
metadata_last[] <- lapply(metadata_last, function(x) ifelse(x == "NA", NA, x))  

# detect numerical columns 
is_numeric_column <- function(x) {
    all(is.na(x) | grepl("^-?\\d+(\\.\\d+)?$", x))
} 
numeric_cols <- sapply(metadata_last, is_numeric_column)
# colnames(metadata.fecal.baseline)[numeric_cols] # the name of detected numerical variables 

# convert numerical column into numerical variable type  
metadata_last[numeric_cols] <- lapply(metadata_last[numeric_cols], as.numeric) 
table(metadata_last$group)

# select groups
metadata_last <- metadata_last %>% filter(group %in% c("HC", "GBA_NMC", "GBA_PD", "sPD"))



###############################
# Control quality 
###############################
# Discard fecal samples based on manual and crocodeel QC validation, performed by Victoria Meslier. 

# Selection QC 
# fecal samples to exclude based on manual and crocodeel QC validation and one sample ASAPh_1114 has failed the sequencing after several attempts
fecal_samples_to_exclude <- c("ASAPh_024","ASAPh_076", "ASAPh_1081", "ASAPh_1082", "ASAPh_1085", "ASAPh_1088", "ASAPh_1097", "ASAPh_1102", "ASAPh_1111", "ASAPh_1113", "ASAPh_1133", "ASAPh_1144", "ASAPh_1148", "ASAPh_1152", "ASAPh_1167", "ASAPh_1168", "ASAPh_229", "ASAPh_647", "ASAPh_650", "ASAPh_661", "ASAPh_664", "ASAPh_681", "ASAPh_811", "ASAPh_818", "ASAPh_821", "ASAPh_822", "ASAPh_835", "ASAPh_870", "ASAPh_878", "ASAPh_969", "ASAPh_972", "ASAPh_973", "ASAPh_982", "ASAPh_983", "ASAPh_990", "ASAPh_642","ASAPh_802", "ASAPh_815", "ASAPh_836", "ASAPh_840", "ASAPh_841", "ASAPh_843", "ASAPh_853", "ASAPh_860", "ASAPh_866", "ASAPh_879", "ASAPh_975", "ASAPh_979", "ASAPh_986", "ASAPh_1107", "ASAPh_1131", "ASAPh_1132", "ASAPh_1141", "ASAPh_1150", "ASAPh_1151", "ASAPh_1114") # 56samples
fecal_samples_to_exclude2 <- c("ASAPh_533", "ASAPh_798", "ASAPh_820", "ASAPh_824", "ASAPh_832", "ASAPh_851")

#creating a variable 'to_exclude_fecal_baseline' that exclude patients for which no stool and/or saliva samples was obtained and samples that do not pass our QC control
metadata_last <- metadata_last %>% 
    mutate(to_exclude_fecal_baseline = case_when(is.na(MGP_ID_fecal_baseline)~"yes", 
                                                 MGP_ID_fecal_baseline %in% fecal_samples_to_exclude ~ "yes", 
                                                 MGP_ID_fecal_baseline %in% fecal_samples_to_exclude2 ~ "yes", 
                                                 TRUE ~ "no"))

metadata.fecal.baseline <- metadata_last %>% 
    filter(to_exclude_fecal_baseline == "no")
rownames(metadata.fecal.baseline) <- metadata.fecal.baseline$MGP_ID_fecal_baseline
dim(metadata.fecal.baseline)



# Variables from clinical codebook (slight differences or variation in name)
selected_var <- c("MGP_ID_fecal_baseline", "cohort", "age", "sex", "partner", "blood_related", "bmi", "group", "fhx", "pd_criteria", "tbi", "appendectomy", "appendectomy_ageatsurgery", "appendectomy_pd_yearsdiff", "coffee", "coffee_amount", "smoking", "smoking_current_duration", "smoking_ex_duration", "smoking_packs", "updrs_i_total", "updrs_ii_total",  "hoehn_yahr",  "updrs_iii_total", "updrs_iv_total", "updrs_total", "moca_raw_total", "education_n_years", "ed_yrs", "moca_total", "scopa_total", "scopa_gastrointestinal", "scopa_urinary", "scopa_cardiovascular", "scopa_thermoregulatory", "scopa_pupillomotor", "scopa_sexual", "hads_depression_total", "hads_anxiety_total", "rbdsq_total", "bdi_total", "upsit_total", "wcs_total", "seborr_dermatitis", "seborr_dermatitis_ageonset", "seborr_dermatitis_pd_datediff", "pd_dbs", "pd_dbs_years", "pd_dbs_target", "pd_lcig", "pd_apomorphine", "pd_onset_calculated", "pd", "pd_dur", "pd_drug_naive", "ledd_total", "diet", "fruit_portion", "vegetables_portion", "milk", "alcohol_intake", "smoking_pd_biome", "physical_activity", "physical_activity_threshold", "out_breath_stairs", "Fruit", "Vegetables", "Oily_Fish", "Fat", "NMES", "DQS_score", "hads_depression_range", "hads_anxiety_range", "bdi_range", "rbd_affected_diff",  "upsit_range", "wcs_threshold", "moca_range", "subthr_park", "subthr_park_threshold", "sleepiness", "sleepiness_prodromal_criteria", "urinary_dysfunction", "urinary_dysfunction_prodromal_criteria", "erectyle_dysfunction", "erectyle_dysfunction_prodromal_criteria", "physical_activity_threshold_prodromal_criteria", "t2_diabetes", "prior_probability", "risk_markers", "prodromal_markers", "clinical_nonmotor_markers", "clinical_motor_markers", "lr_threshold", "estimated_lr", "estimated_probability", "prob_prodromalPD") 
metadata.fecal.baseline <- metadata.fecal.baseline %>% 
    select(all_of(selected_var))



###############################
# Merging fecal MSP matrix (gut & oral) at 18M reads per sample for further analysis 
###############################
# loading gut MSP matrix (output of meteor2)
msp.gut.531ind = read.delim("./data/ASAP_humain_vs_10_4_igc2_smart_shared_reads_531ind_fecal_final_14122023.downHQ.18000000.fpkm.mgs.100.tsv",header = T, stringsAsFactors = F, row.names = 1)
dim(msp.gut.531ind) #1990 MSP and 531 samples
msp.gut.464ind = msp.gut.531ind[,match(rownames(metadata.fecal.baseline), colnames(msp.gut.531ind))] # 464

# loading oral MSP matrix (output of meteor2)
msp.oral.531ind = read.delim("./data/ASAP_humain_vs_8_4_oral_smart_shared_reads_531ind_fecal_final_15122023.downHQ.18000000.fpkm.mgs.100.tsv",header = T, stringsAsFactors = F, row.names = 1)
dim(msp.oral.531ind) #853 MSP and 531 samples
msp.oral.464ind = msp.oral.531ind[,match(rownames(metadata.fecal.baseline), colnames(msp.oral.531ind))] # 464

# checking colnames
all(colnames(msp.oral.464ind) == colnames(msp.gut.464ind)) # T

# loading species reference (taxo), file build by Emmanuelle Le Chatelier & Florian Plaza Oñate
MSP_set = get(load("./data/MSP_set_ref_oss_gss_status_612_ref_paper_20220622.RData")) 

# set origins of the reference MSP signal (oral or gut)
MSP_set_gut <- MSP_set[MSP_set$ref_MSP == "gut",] # 1947
MSP_set_oral <- MSP_set[MSP_set$ref_MSP == "oral",] # 794

# merge oral and gut species tables to get unique MGS signal table
msp_fecal_baseline_464ind <- rbind(msp.gut.464ind[rownames(msp.gut.464ind) %in% rownames(MSP_set_gut),],
                                   msp.oral.464ind[rownames(msp.oral.464ind) %in% rownames(MSP_set_oral),]) 
dim(msp_fecal_baseline_464ind) # 2741 MSP 464 samples

# Filter MSP with null abundance 
msp_fecal_baseline_464ind <- msp_fecal_baseline_464ind[rowSums(msp_fecal_baseline_464ind) >0, colSums(msp_fecal_baseline_464ind) > 0] # 2054 MSP left and 464ind


# add msp_1786 in msp abundance matrix, with abundance = 0 for each individual in order to avoid warning when performing potential functional analysis (msp_1786 was present in completude but was absent in abundance matrix)
msp_1786 <- rep(0, ncol(msp_fecal_baseline_464ind))
msp_fecal_baseline_464ind <- rbind(msp_fecal_baseline_464ind, msp_1786) 
rownames(msp_fecal_baseline_464ind)[nrow(msp_fecal_baseline_464ind)] <- "msp_1786"  

dim(msp_fecal_baseline_464ind) # 2055 MSP and 464ind 


###############################
# Making matrices tables for higher taxonomical ranks (not used)
###############################
#computing norm composition
msp_fecal_baseline_464ind_norm = msp_fecal_baseline_464ind
for(i in 1:ncol(msp_fecal_baseline_464ind_norm)){
    msp_fecal_baseline_464ind_norm[,i] <- msp_fecal_baseline_464ind_norm[,i]/colSums(msp_fecal_baseline_464ind_norm)[i]}
# colSums(msp_fecal_baseline_464ind_norm) #ok checked



###############################
# alpha diversity indices determination
###############################
#adding MSP richness from merged data
all(colnames(msp_fecal_baseline_464ind) == rownames(metadata.fecal.baseline)) #T
metadata.fecal.baseline$MSP_richness_fecal_baseline = colSums(msp_fecal_baseline_464ind>0)

# adding GC from gut catalog
GC_richness = read.delim("./data/ASAP_humain_vs_10_4_igc2_smart_shared_reads_531ind_fecal_final_14122023.richness.txt", sep=" ", row.names=1, stringsAsFactors = F) # contains GC richness computed for each sample

metadata.fecal.baseline$GC_richness_gutcat_fecal_baseline = GC_richness[match(metadata.fecal.baseline$MGP_ID_fecal_baseline, rownames(GC_richness)), "down_18000000_una"]

#adding alpha div indices
alpha.div <- as.data.frame(diversity(t(msp_fecal_baseline_464ind),index = "shannon"))
colnames(alpha.div)[1] <- "shannon"
alpha.div$simpson <- diversity(t(msp_fecal_baseline_464ind),index = "simpson")

data <- apply(t(msp_fecal_baseline_464ind)>0, 1, sum)
alpha.div$evenness <- diversity(t(msp_fecal_baseline_464ind), index = "simpson")/log(data)

# adding to metadata
metadata.fecal.baseline$shannon_fecal_baseline <- alpha.div[match(metadata.fecal.baseline$MGP_ID_fecal_baseline, rownames(alpha.div)), "shannon"]
metadata.fecal.baseline$simpson_fecal_baseline <- alpha.div[match(metadata.fecal.baseline$MGP_ID_fecal_baseline, rownames(alpha.div)), "simpson"]
metadata.fecal.baseline$evenness_fecal_baseline <- alpha.div[match(metadata.fecal.baseline$MGP_ID_fecal_baseline, rownames(alpha.div)), "evenness"]


###############################
# Functional data preparation - no core option - code from Florence Thirion
###############################
# Loading modules definitions (file available at https://github.com/metagenopolis/meteor/releases/tag/2.0.18, under the name of : module.feather)
load(file = "./data/all_modules_definition_GMM_GBM_KEGG_107_jan2024.RData") 
colnames(all_modules_definition_GMM_GBM_KEGG_107)[7] <- "GBM_function"
all_modules_definition_GMM_GBM_KEGG_107 <- as.data.frame(all_modules_definition_GMM_GBM_KEGG_107)

# Load completude (outputs of meteor2) and fusion completude 
completude.gut = fread("./data/completude_ind_mgs_mod_data_downHQ_18000000fpkm_ftmt.txt", header = FALSE) 
completude.gut = as.data.frame(completude.gut[,c(1:4)])
colnames(completude.gut) = c("ind", "mod", "msp", "c0")

completude.oral = fread("./data/completude_ind_mgs_mod_data_downHQ_18000000fpkm_ftmt.txt", header = FALSE)
completude.oral = as.data.frame(completude.oral[,c(1:4)])
colnames(completude.oral) = c("ind", "mod", "msp", "c0")

# merge oral and gut functional data and keeping only modules with completion > 0.9
completude.nocore = rbind(subset(completude.gut, c0 >= 0.9),
                          subset(completude.oral, c0 >= 0.9))
completude.nocore$unique = paste0(completude.nocore$ind, completude.nocore$mod, completude.nocore$msp)
completude.nocore = completude.nocore[!duplicated(completude.nocore$unique),]

# filtering to keep only fecal baseline ind
completude.nocore.fecal.baseline <- completude.nocore %>% filter(ind %in% metadata.fecal.baseline$MGP_ID_fecal_baseline)
completude.nocore.fecal.baseline$ind %>% unique() %>% length() #464individuals OK


# making the functional abundance data
ftmt.nocore.fecal.baseline = matrix(0, ncol = length(unique(completude.nocore.fecal.baseline$ind)), nrow = length(unique(completude.nocore.fecal.baseline$mod)))
rownames(ftmt.nocore.fecal.baseline) = unique(completude.nocore.fecal.baseline$mod)
colnames(ftmt.nocore.fecal.baseline) = unique(completude.nocore.fecal.baseline$ind)
ftmt.nocore.fecal.baseline = as.data.frame(ftmt.nocore.fecal.baseline)

for (i in 1:length(unique(completude.nocore.fecal.baseline$ind))){
    my_ind = unique(completude.nocore.fecal.baseline$ind)[i]
    my_subset = subset(completude.nocore.fecal.baseline, ind == my_ind) # recover all lines (=mod*msp) for this individual
    my_subset = split(my_subset, my_subset$mod) # make a list for each module containing msp that carry this module for this individual
    tmp_res = sapply(my_subset, function(x) sum(msp_fecal_baseline_464ind[x[["msp"]], my_ind])) # module = sum of abundance of MSP carrying that module 
    ftmt.nocore.fecal.baseline[names(tmp_res), my_ind] = tmp_res # report results in final table
}

dim(ftmt.nocore.fecal.baseline) #357 modules and 464 individuals

# checking : now no NA 
# ftmt.nocore.fecal.baseline["MGB031", "ASAPh_576"]
# ftmt.nocore.fecal.baseline["M00307", "ASAPh_576"]
# ftmt.nocore.fecal.baseline["MF0061", "ASAPh_576"]




# There are variables that are not well recognize as numerical variables, here we convert them into numerical type. 
# Function to detect numerical columns with decimals (taking into account missing values) 
is_numeric_column <- function(x) {
    all(is.na(x) | grepl("^-?\\d+(\\.\\d+)?$", x))
} 

numeric_cols <- sapply(metadata.fecal.baseline, is_numeric_column)
# colnames(metadata.fecal.baseline)[numeric_cols] # the name of detected numerical variables 

# again, convert numerical column into numerical variable type  
metadata.fecal.baseline[numeric_cols] <- lapply(metadata.fecal.baseline[numeric_cols], as.numeric)
# str(metadata.fecal.baseline)


###############################
# Save RData 
###############################
save(metadata.fecal.baseline, file = "./data/ASAPh_metadata_fecal_baseline_464ind_10102025_victoria.RData")
save(msp_fecal_baseline_464ind, MSP_set, completude.nocore.fecal.baseline, ftmt.nocore.fecal.baseline, all_modules_definition_GMM_GBM_KEGG_107, file = "./data/ASAPh_microbial_data_fecal_baseline_464ind_10102025_victoria.RData")

