# Data management and preparation for the fecal baseline analysis with 464 individuals 

###############################
# Merging fecal MSP matrix (gut & oral) at 18M reads per sample for further analysis 
###############################
# loading gut MSP matrix 
msp.gut.531ind = read.delim("./data/ASAP_humain_vs_10_4_igc2_smart_shared_reads_531ind_fecal_final_14122023.downHQ.18000000.fpkm.mgs.100.tsv",header = T, stringsAsFactors = F, row.names = 1)
dim(msp.gut.531ind) #1990 MSP and 531 samples
msp.gut.464ind = msp.gut.531ind[,match(rownames(metadata.fecal.baseline), colnames(msp.gut.531ind))] # 464

# loading oral MSP matrix 
msp.oral.531ind = read.delim("./data/ASAP_humain_vs_8_4_oral_smart_shared_reads_531ind_fecal_final_15122023.downHQ.18000000.fpkm.mgs.100.tsv",header = T, stringsAsFactors = F, row.names = 1)
dim(msp.oral.531ind) #853 MSP and 531 samples
msp.oral.464ind = msp.oral.531ind[,match(rownames(metadata.fecal.baseline), colnames(msp.oral.531ind))] # 464

# checking colnames
all(colnames(msp.oral.464ind) == colnames(msp.gut.464ind)) # T

# loading species reference (taxo) 
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


# add msp_1786 in msp abundance matrix, with abundance = 0 for each individual in order to avoid warning when performing potential functional analysis
msp_1786 <- rep(0, ncol(msp_fecal_baseline_464ind))
msp_fecal_baseline_464ind <- rbind(msp_fecal_baseline_464ind, msp_1786) 
rownames(msp_fecal_baseline_464ind)[nrow(msp_fecal_baseline_464ind)] <- "msp_1786"  

dim(msp_fecal_baseline_464ind) # 2055 MSP and 464ind 