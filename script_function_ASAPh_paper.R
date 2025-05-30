

# function that compute the regression from a dataframe 
reg_plot <- function(df_melt, x, y, title = NULL, xlab = NULL, ylab = NULL, color_var = NULL, point_shape = 19, formula.pos = 0.9){
    formula <- y ~ x
    p <- ggplot(data = df_melt, aes(x = !!as.symbol(x), y = !!as.symbol(y)))
    
    if (color_var %in% colnames(df_melt)){
        p <- p + geom_point(aes(color = !!as.symbol(color_var)), size = 4, shape = point_shape) +
            scale_color_manual(values = c("black", "red"))
    } 
    else{ 
        p <- p + geom_point(color = color_var, size = 4, shape = point_shape)
    }
    p <- p + 
        theme_classic() +
        geom_smooth(method = "lm", se = F, formula = formula, linetype = 3, linewidth=1, color = "black") +
        stat_poly_eq(use_label(c("eq", "R2")), label.y = formula.pos, size = 7, output.type = "expression") +

        xlab(xlab) +
        ylab(ylab) + 
        theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 20))
    return(p)
}

# function create a dataframe with 2 columns : group and cliffdelta (absolute value) 
get_cliff_delta <- function(df_merge_compa, cat){
    df_merge_compa <- df_merge_compa %>% filter(type == cat)
    
    cliff <- c(df_merge_compa$HC.GBA_NMC_Cliffdelta, df_merge_compa$HC.PD_Cliffdelta) %>% abs()
    
    if (cat == "CE"){
        cliff_status <- c(df_merge_compa$HC.GBA_NMC_status_Cliffdeltacorrected,df_merge_compa$HC.PD_status_Cliffdeltacorrected)
    } 
    if (cat == "CD"){
        cliff_status <- c(df_merge_compa$HC.GBA_NMC_status_Cliff_depleted, df_merge_compa$HC.PD_status_Cliff_depleted)
    }
    if (cat == "NE"){
        cliff_status <- c(df_merge_compa$HC.GBA_NMC_status_Cliff_depleted, df_merge_compa$HC.PD_status_Cliffdeltacorrected)
    }    
    if (cat == "ND"){
        cliff_status <- c(df_merge_compa$HC.GBA_NMC_status_Cliffdeltacorrected, df_merge_compa$HC.PD_status_Cliff_depleted)
    }
    df_final <- data.frame(cliff, cliff_status)
    return(df_final)
}

# function that performs Student test between 2 groups of cliffdelta values
compute_cliff_ttest <- function(df_cliff){ 
    res_test <- data.frame()
    formula <- as.formula(paste("cliff", "cliff_status", sep = "~"))
    mean_cliff <- aggregate(formula, data = df_cliff, FUN = mean)
    direction <- ifelse(mean_cliff[mean_cliff$cliff_status == "GBA_NMC", "cliff"] < mean_cliff[mean_cliff$cliff_status == "PD", "cliff"], "less", "greater") 
    print(direction)
    
    res_test <- t_test(data = df_cliff, formula = formula, alternative = direction, var.equal = T)
    
    colnames(res_test)[1] <- "y"
    return(res_test)
}


# function that plots a formated table containing the pvalues of Student test 
get_tab_pval_ttest <- function(list_CE_prop, list_CD_prop){
    my_compa <- combn(vec_group, m = 2)
    table_pval <- data.frame(matrix(nrow = 2, ncol = ncol(my_compa)))
    rownames(table_pval) <- c("PD enriched", "PD depleted")
    
    for (i in 1:ncol(my_compa)){
        grp1 = my_compa[1,i]
        grp2 = my_compa[2,i]
        # for checking : 
        colnames(table_pval)[i] <- paste0(grp1,"\n", "vs\n", grp2)
        
        # Enriched
        uni_CE <- ifelse(mean(list_CE_prop[[grp1]]) < mean(list_CE_prop[[grp2]]), "less", "greater")
        print(uni_CE)
        res_student_CE <- t.test(as.numeric(list_CE_prop[[grp1]]), as.numeric(list_CE_prop[[grp2]]), alternative=uni_CE, var.equal=T)
        table_pval[1, i] <- signif(res_student_CE$p.value, 2) 
        
        # Depleted
        uni_CD <- ifelse(mean(list_CD_prop[[grp1]]) < mean(list_CD_prop[[grp2]]), "less", "greater") 
        print(uni_CD)
        res_student_CD <- t.test(as.numeric(list_CD_prop[[grp1]]), as.numeric(list_CD_prop[[grp2]]), alternative=uni_CD, var.equal=T)
        table_pval[2, i] <- signif(res_student_CD$p.value, 2)
    }
    
    # format pvalues into scientific notation  
    table_pval <- table_pval %>% dplyr::mutate_if(is.numeric, list(~format(., scientific = T)))
    return(table_pval)
} 

# function that computes proportion of subtype specie relative to the group 
compute_prop <- function(vec_group, metadata, abund_mat, vec_msp){ 
    grp_sub_prop <- list()
    for (grp in vec_group){ 
        ind <- metadata %>% dplyr::filter(group_bis == grp) %>% pull(MGP_ID_fecal_baseline)
        # count nb of detected MSP in samples of each group  
        grp_abund <- abund_mat[, ind]
        grp_count <-  apply(grp_abund, 2, function(x) {length(which(x>0))})
        
        # count nb of specific MSP in samples of each group, and compute its proportion
        grp_sub_abund <- abund_mat[vec_msp, ind]
        grp_sub_count <- apply(grp_sub_abund, 2, function(x) {length(which(x>0))})
        
        
        prop <- grp_sub_count/grp_count*100
        grp_sub_prop[[grp]] <- prop
    }
    return(grp_sub_prop)
}

# function that computes abundance of subtype specie relative to the group 
compute_abund <- function(vec_group, metadata, abund_mat, vec_msp){ 
    grp_sub_abund <- list()
    for (grp in vec_group){
        ind <- metadata %>% dplyr::filter(group_bis == grp) %>% pull(MGP_ID_fecal_baseline)
        # sum MSP abundance in samples of each group 
        grp_abund <- abund_mat[, ind]
        grp_sum_abund <- apply(grp_abund, 2, sum)
        
        # sum specific MSP abundance in samples of each group and compute its abundance in samples relative to the group
        sub_abund <-  abund_mat[vec_msp, ind]
        sub_sum_abund <- apply(sub_abund, 2, sum)
        
        rel_abund <- sub_sum_abund/grp_sum_abund*100
        grp_sub_abund[[grp]] <- rel_abund
    } 
    return(grp_sub_abund)
}  


# function that construct barplot
show_barplot <- function(df_melt, x, y, fill, my_color, ylab = NULL, title = NULL, limits = NULL, pas = 1, show_legend){
    p <- ggbarplot(df_melt, x = x, y = y, fill = fill, position = position_dodge(), ylab = ylab, title = title) + 
        scale_fill_manual(values = my_color) + 
        theme_classic()+
        theme(legend.title=element_blank(), 
              legend.position = show_legend, 
              plot.title = element_text(size=18),
              axis.title.x = element_blank(),
              axis.title = element_text(size = 18, color = "black"),
              axis.text = element_text(size = 18, color = "black"), 
              plot.margin=unit(c(0.1,0,0,0.1), "cm")) 
    if (!is.null(limits)){
        p <- p + expand_limits(y = limits) + scale_y_continuous(breaks = seq(limits[1], limits[2], pas))  
    }
    return(p)
}


show_lm_eq <- function(data, x, y){
    formula <- as.formula(paste(y, x, sep = " ~ "))
    m <- lm(formula = formula, data)
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                     list(a = format(coef(m)[[1]], digits = 3), 
                          b = format(coef(m)[[2]], digits = 3), 
                          r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));                 
}


show_tab <- function(tab){ 
    final_tab <- ggpubr::ggtexttable(tab, theme = ttheme(base_style  = "classic", base_size = 17)) %>% 
        tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>% 
        tab_add_hline(at.row = 3, row.side = "bottom", linewidth = 3, linetype = 1) + theme(plot.margin = unit(c(0,0,0,0), 'lines'))
    return(final_tab)
}


# function that compute the mean of proportion/abundance of enriched/depleted MSP in each group 
show_mean_tab <- function(list_grp1, list_grp2, vec_cat){
    # compute mean for each sublist (species) from the list (group)
    mean_grp1 <- lapply(list_grp1, mean) 
    mean_grp2 <- lapply(list_grp2, mean)
    
    # create tab of mean values : concatenate the rows from the 2 df
    tab_mean <- do.call(rbind.data.frame, list(mean_grp1, mean_grp2)) 
    tab_mean$cat <- vec_cat
    tab_mean_melt <- reshape2::melt(tab_mean)
    tab_mean_melt <- tab_mean_melt %>% mutate(cat = factor(cat, levels = vec_cat))
    return(tab_mean_melt)
}


# function that computes proportion of subtype specie relative to the group 
compute_prop_subtype <- function(vec_sub_group, specific_ind, abund_mat){ 
    list_sub_prop <- list()
    for (sub in vec_sub_group){ 
        # count nb of detected MSP in samples of each group  
        grp_abund <- abund_mat[, specific_ind]
        grp_count <-  apply(grp_abund, 2, function(x) {length(which(x>0))})
        
        # count nb of specific MSP in samples of each group, and compute its proportion
        sub_abund <- abund_mat[get(sub), specific_ind]
        sub_count <- apply(sub_abund, 2, function(x) {length(which(x>0))})
        
        check <- all(names(sub_count) == names(grp_count))
        print(check)
        prop <- sub_count/grp_count*100
        list_sub_prop[[sub]] <- prop
    }
    return(list_sub_prop)
}

# function that computes abundance of subtype specie relative to the group
compute_abund_subtype <- function(vec_sub_group, specific_ind, abund_mat){ 
    list_sub_abund <- list()
    for (sub in vec_sub_group){
        # sum MSP abundance in samples of each group 
        grp_abund <- abund_mat[, specific_ind]
        grp_sum_abund <- apply(grp_abund, 2, sum)
        
        # sum specific MSP abundance in samples of each group and compute its abundance in samples relative to the group
        sub_abund <-  abund_mat[get(sub), specific_ind]
        sub_sum_abund <- apply(sub_abund, 2, sum)
        
        check <- all(names(sub_sum_abund) == names(grp_sum_abund)) 
        print(check)
        rel_abund <- sub_sum_abund/grp_sum_abund*100
        list_sub_abund[[sub]] <- rel_abund
    } 
    return(list_sub_abund)
}   

# function that compute Student test
compute_ttest <- function(list1, list2, grp, grp_label){ 
    vec_sub_group <- c("CE", "NE", "CD", "ND")
    df1 <- as.data.frame(list1)
    df1[, grp] <- grp_label[1]
    df2 <- as.data.frame(list2)
    df2[, grp] <- grp_label[2]
    
    df <- rbind(df1, df2)
    df[, grp] <- factor(df[, grp], levels = grp_label)
    print(df)
    
    res_test <- data.frame()
    for (sub in vec_sub_group){ 
        formula <- as.formula(paste(sub, grp, sep = " ~ "))
        direction <- ifelse(mean(list1[[sub]]) < mean(list2[[sub]]), "less", "greater")
        print(mean(list1[[sub]]))
        print(mean(list2[[sub]]))
        #print(formula)
        print(paste(sub, grp, sep = " ~ "))
        print(direction)
        
        res <- t_test(data = df, formula = formula, alternative = direction, var.equal = T)
        res_test <- res_test %>% rbind(res)
    }
    colnames(res_test)[1] <- "y"
    return(res_test)
}

# FUNCTIONS QUARTILES
# function that sort proportion/abundance in ascending way for enriched subtype specie
# in descending way for depleted subtype specie and split into quartile for PD
get_quartile_PD <- function(data_list, msp_subgroup = "E"){
    # ascending sort of PD individuals according to the proportion/abundance of enriched MSP 
    if (msp_subgroup == "E"){
        sorted <- sort(data_list)
    }
    if (msp_subgroup == "D"){
        sorted <- sort(data_list, decreasing = T)
    }
    Q1 <- sorted[1:68]
    Q2 <- sorted[69:136]
    Q3 <- sorted[137:203]
    Q4 <- sorted[204:271]
    
    quartile <- list(Q1, Q2, Q3, Q4)
    names(quartile) <- c("Q1", "Q2", "Q3", "Q4")
    return(quartile)
}



get_quartile_NMC2 <- function(data_list, vec_msp){
  n <- length(data_list)
  print(n)
  
  if (vec_msp %in% c("CE", "NE")){ 
    sorted <- sort(data_list)
    # if n is pair 
    if (n %% 2 == 0){ 
      index_median1 <- n / 2
      index_median2 <- n / 2 + 1
      below <- sorted[1:index_median1]
      above <- sorted[index_median2:n]
      # if n is impair 
    } else {
      index_median <- (n + 1) / 2
      # the median belong to the first group 
      below <- sorted[1:(index_median-1)]
      above <- sorted[index_median:n]
    }
  } else if (vec_msp %in% c("CD", "ND")){
    sorted <- sort(data_list, decreasing = T)
    # if n is pair 
    if (n %% 2 == 0){ 
      index_median1 <- n / 2
      index_median2 <- n / 2 + 1
      below <- sorted[1:index_median1]
      above <- sorted[index_median2:n]
      # if n is impair 
    } else {
      index_median <- (n + 1) / 2
      print(index_median)
      # the median belong to the first group 
      below <- sorted[1:index_median]
      above <- sorted[(index_median+1):n]
    }
  }
  return(list("below" = below, "above" = above))
}  


# function that sort proportion/abundance in ascending way for enriched subtype specie
# in descending way for depleted subtype specie and split into quartile for HC
get_quartile_HC <- function(data_list, msp_subgroup = "E"){
    # ascending sort of HC individuals according to the proportion/abundance of enriched MSP 
    if (msp_subgroup == "E"){
        sorted <- sort(data_list)
    }
    if (msp_subgroup == "D"){
        sorted <- sort(data_list, decreasing = T)
    }
    Q1 <- sorted[1:38]
    Q2 <- sorted[39:75]
    Q3 <- sorted[76:112]
    Q4 <- sorted[113:150]
    
    quartile <- list(Q1, Q2, Q3, Q4)
    names(quartile) <- c("Q1", "Q2", "Q3", "Q4")
    return(quartile)
}



# function that compute mean using a list 
compute_mean <- function(list){ 
    mean_list <- lapply(list, mean)
    df <- mean_list %>% as.data.frame() %>% reshape2::melt() 
    return(df)
}



# function that take the matrices from rcorr output and arrange 
get_corr_and_pval <- function(rcorr_out, abund_mat_grp){ 
    # get the matrix of correlation and its pvalues 
    mat_cor <- rcorr_out$r
    mat_corP <- rcorr_out$P
    
    # prevent redundancy in matrices :  keep MSP in rows and clinical variables in columns
    sub_mat_cor <- mat_cor[1:nrow(abund_mat_grp), (nrow(abund_mat_grp)+1):ncol(mat_cor)]
    sub_mat_corP <- mat_corP[1:nrow(abund_mat_grp), (nrow(abund_mat_grp)+1):ncol(mat_corP)] 
    
    # reshape matrices correlation and pvalue matrices
    sub_mat_cor_melt <- reshape2::melt(sub_mat_cor)
    colnames(sub_mat_cor_melt) <- c("msp_id", "variable", "corr")
    sub_mat_corP_melt <- reshape2::melt(sub_mat_corP)
    colnames(sub_mat_corP_melt) <- c("msp_id", "variable", "p_corr") 
    
    # merge correlation and pvalue matrices
    merge_mat_cor <- merge(sub_mat_cor_melt, sub_mat_corP_melt, by = c("msp_id", "variable"))
    
    # define subgroup (CE, NE, CD, ND)
    merge_mat_cor <- merge_mat_cor %>% 
        mutate(type = case_when(msp_id %in% CE ~ "CE", 
                                msp_id %in% NE ~ "NE", 
                                msp_id %in% CD ~ "CD", 
                                msp_id %in% ND ~ "ND"))
    
    # MSP that belong to CE and NE subgroups are defined as PD enriched, otherwise defined as PD depleted
    merge_mat_cor <- merge_mat_cor %>% 
        mutate(cat = ifelse(type %in% c("CE", "NE"), "PD_enriched", "PD_depleted")) %>%
        mutate(cat = factor(cat, levels = c("PD_enriched", "PD_depleted")))
    
    # create new variable that indicates if correlation is positive or negative
    merge_mat_cor <- merge_mat_cor %>% mutate(corr_dir = ifelse(corr>0, "pos", "neg")) #%>% 
    #mutate(corr_dir = factor(corr_dir, levels = c("pos", "neg")))
    
    return(merge_mat_cor)
}



# function that compute mean and standard deviation using two subset of individuals  
compute_mean_sd <- function(subset1, subset2, metadata, name_grp){
  # get the metadata for selected individuals and compute its mean and sd for each variable
  subset1_data <- metadata %>% filter(rownames(metadata) %in% names(subset1))
  subset1_mean <- apply(subset1_data, 2, function(column) round(mean(column, na.rm = T), 1))
  subset1_sd <- apply(subset1_data, 2, function(column) round(sd(column), 1))
  # for checking purpose 
  subset1_with_na <- colnames(subset1_data)[apply(subset1_data, 2, function(column) any(is.na(column)))]
  
  subset2_data <- metadata %>% filter(rownames(metadata) %in% names(subset2))
  subset2_mean <- apply(subset2_data, 2, function(column) round(mean(column, na.rm = T), 1))
  subset2_sd <- apply(subset2_data, 2, function(column) round(sd(column, na.rm = T), 1))
  subset2_with_na <- colnames(subset2_data)[apply(subset2_data, 2, function(column) any(is.na(column)))]
  
  # print("Subset1, NA variables : ", subset1_with_na)
  # print("Subset2, NA variables : ", subset2_with_na)
  
  subset1_name <- paste0(name_grp, "_L")
  subset2_name <- paste0(name_grp, "_M")  
  list_mean <- setNames(list(subset1_mean, subset2_mean), 
                        c(subset1_name, subset2_name))
  
  list_sd <- setNames(list(subset1_sd, subset2_sd), 
                      c(subset1_name, subset2_name))
  
  
  return(list("mean" = list_mean, "sd" = list_sd,  "data1" = subset1_data, "data2" = subset2_data)) 
}


# function that perform Student test using t_test function and 
# adjust the symbol of significance to pvalues (different from the those of t_test)
compute_pval <- function(df){
    res_test <- data.frame()
    for (i in 1:(ncol(df)-1)){ 
        variable <- colnames(df)[i] 
        formula <- as.formula(paste(variable, "quartile", sep = " ~ ")) 
        mean1 <- mean(df[which(df$quartile == "Q1") ,variable], na.rm = T)
        mean2 <- mean(df[which(df$quartile == "Q4") ,variable], na.rm = T)
        direction <- ifelse(mean1 < mean2, "less", "greater")
        print(direction)
        
        res <- t_test(data = df, formula = formula, alternative = direction, var.equal = T) 
        res <- res %>% mutate(p.signif = case_when(p < 0.0001 ~ "****",
                                                   p < 0.001 ~ "***",
                                                   p < 0.01 ~ "**",
                                                   p < 0.05 ~ "*", 
                                                   p >= 0.05 ~ ""))
        res_test <- res_test %>% rbind(res)
    }
    colnames(res_test)[1] <- "y"
    return(res_test)
  
}



################################################################################################
#Function used to produce Figure 4#
################################################################################################
# Firstly a function to generate correlation table.
cor_176_clinic <- function(list_df_cor, groups, data, metadata, clinic_var, name_var) {
  
  # Looping over groups 
  for (gp in groups) {
    
    names_nw_df <- paste0("correlation_", gp)
    
    intermediare_df <- data.frame()
    
    ## keeping the current group during each iteration.
    data_filter <- data %>%
      filter(group == gp) %>%
      select(-group)
    
    # Looping over clinical variable
    for (clin_var in clinic_var) {
      
      ## from data filter from groups keep only the current clinical variable
      data_filter_clinic <- data_filter %>%
        merge(metadata[clin_var, drop = FALSE], by = "row.names") %>%
        column_to_rownames(var = "Row.names")
      
      ### setdiff is for only keep the msp name and remove the clinical var added previously
      for (col in setdiff(names(data_filter_clinic), clin_var)) {
        
        ### if there is as much NA than there is row set the correlation and pvalue to NA
        if (sum(is.na(data_filter_clinic[[clin_var]])) == nrow(data_filter)){
          
          intermediare_df <- rbind(intermediare_df, as.data.frame(setNames(
            list(col, clin_var, "NA", "NA"),
            c(name_var, "clinical_variable", "correlation", "pvalue"))))
          
          ### otherwise calculation of the spearman correlation value
        } else {
          
          cor_res <- cor.test(data_filter_clinic[[col]], data_filter_clinic[[clin_var]], method = "spearman")
          
          intermediare_df <- rbind(intermediare_df, as.data.frame(setNames(
            list(col, clin_var, cor_res$estimate, cor_res$p.value),
            c(name_var, "clinical_variable", "correlation", "pvalue")
          )))
          
        }
        
      }
    }
    
    ### store the intermediaire dataframe to a list with corresponding name 
    list_df_cor[[names_nw_df]] <- intermediare_df  
  }
  
  return(list_df_cor)  
}


# Then, a function for generating tables in a format easier to use.

### to use in order to transform initial dataframe into new format needed to plot.
### format expected is : str(PD)
# 'data.frame':	4844 obs. of  12 variables:
#   $ msp_id           : chr  "msp_0005" "msp_0013" "msp_0014" "msp_0015" ...
# $ msp_name         : chr  "msp_0005_Escherichia coli" "msp_0013_Ruminococcus bicirculans" "msp_0014_Eisenbergiella tayi" "msp_0015_Agathobacter faecis" ...
# $ clinical_variable: chr  "age" "age" "age" "age" ...
# $ correlation      : num  0.12478 -0.00312 0.16579 -0.03028 -0.00917 ...
# $ pvalue           : num  0.0401 0.95929 0.00623 0.61972 0.88051 ...
# $ ref_MSP          : chr  "gut" "gut" "gut" "gut" ...
# $ subgroup         : chr  "CE" "NE" "CE" "NE" ...
# $ prev_PD          : num  234 150 197 200 189 139 174 140 179 179 ...
# $ genus            : chr  "Escherichia" "Ruminococcus 2" "Eisenbergiella" "Agathobacter" ...
# $ family           : chr  "Enterobacteriaceae" "Ruminococcaceae" "Lachnospiraceae G" "Lachnospiraceae C" ...
# $ phylum           : chr  "Proteobacteria" "Firmicutes" "Firmicutes" "Firmicutes" ...
# $ prev_total       : num  390 290 318 360 352 212 321 225 264 301
#### mandatory columns are msp_id, clinical_variable, correlation, pvalue, subgroup

prep_data_fig4 <- function(names_dataframe, # c() c("PD", "GBA_NMC", "HC"),
                           numeric_columns, # correlation, pvalue
                           names_dataframe_filter, # c() c("PD_filter", "GBA_NMC_filter", "HC_filter"),
                           condition_ifelse, # ifelse(data$subgroup %in% c("CE", "ND"), "enriched", "depleted"),
                           clin_var_select, # c() clin_var_select <- c("age", "bdi_total", "bmi","DQS_score", "hads_anxiety_total", "hads_depression_total",
                           # "hoehn_yahr", "moca_total", "rbdsq_total", "scopa_cardiovascular", "scopa_gastrointestinal", "scopa_pupillomotor", "scopa_sexual", "scopa_thermoregulatory", "scopa_total", "scopa_urinary", "updrs_i_total", "updrs_ii_total", "updrs_iii_total", "updrs_total", "upsit_total", "wcs_total", "pd_dur",
                           # "ledd_total"),
                           
                           creating_df, # as.data.frame(matrix(ncol = 4))
                           colnames_df_final #c("clinical_variable", "type", "correlation_pos", "correlation_neg")
                           
){
  for (data_ in names_dataframe){
    data <- get(data_)
    
    name <- paste(data_, "filter", sep = "_")
    
    assign(name, data %>%
             mutate_at(vars(all_of(numeric_columns)), as.numeric) %>%
             filter(pvalue < 0.05))
  }
  
  for(data_ in names_dataframe_filter) {
    data <- get(data_)
    data$type <- eval(condition_ifelse)
    assign(data_, data)
  }
  
  final_df <- setNames(lapply(names_dataframe_filter, function(x) creating_df), names_dataframe_filter)
  
  final_df <- Map(function(df){colnames(df) <- colnames_df_final; df}, final_df)
  
  for(data_ in names_dataframe_filter){
    data <- get(data_)
    for (clin_var in clin_var_select){
      
      
      reshaped_data <- data %>%
        filter(clinical_variable == clin_var) %>%
        group_by(clinical_variable, type) %>%
        summarise(
          correlation_pos = sum(correlation > 0),
          correlation_neg = sum(correlation < 0)) %>%
        distinct(.)
      
      final_df[[data_]] <- rbind(final_df[[data_]], reshaped_data)
      
    }
  }
  
  final_df <- lapply(final_df, function(df) df[-1, ])
  
  final_df <- lapply(final_df, function(df){
    df$chi2 <- sapply(1:nrow(df), function(i){
      chisq.test(c(df$correlation_pos[i], df$correlation_neg[i]), p = c(0.5, 0.5))$p.value
    })
    df$chi2 <- formatC(df$chi2, format = "e", digits = 2)
    df$chi2 <- as.numeric(df$chi2)
    return(df)
  })
  
  return(final_df)
}



### Function that would be use in the plot function.

tobind_row <- function(tobind, rows) {
  for (i in seq_along(rows)) {
    tobind[i, ] <- rows[[i]]
  }
  return(tobind)
}

# Then a function to create the Figure 4
graph_figure4_paper <- function(names_df_tobind, # c() 
                                format_dataframe, #as.data.frame(matrix(ncol = x, nrow = y))
                                rows_tobind, # list() 
                                numeric_columns, # c() 
                                list_clinical_variable = NULL, # c() 
                                final_df, # dataframe
                                clinical_var_col, # string
                                to_exclude, # c() 
                                pos_val, # string
                                neg_val, # string
                                asterisk_pos, #  list()
                                asterisk_neg, # list()
                                name_type_cor_col, # string
                                name_type_cor_val, # string
                                rename_pos_val, # string
                                rename_neg_val, # string
                                ifelse_condition_y, # ifelse()
                                colors_fill, #c()
                                keep_in_legend, # c()
                                labels_legend, # c()
                                which_kind_plot, # c()
                                number_of, # string
                                label_correlation, # c()
                                limits_by, # list()
                                limits_, # list()
                                y_legend, # list()
                                labs_plot, # list("PD_filter" = labs(x = "", y = y_name, fill = ""),"HC_filter" = labs(x = "Clinical variables", y = y_name, fill = ""), "GBA_NMC_filter" = labs(x = "", y = y_name, fill = ""))
                                theme_size, # list()
                                size_asterisk, # numeric
                                angle_asterisk, # numeric
                                x_legend, # numeric
                                size_legend, # numeric
                                letter_fig = NULL, # list()
                                theme_plot_title, #list()
                                theme_health_disease_1,
                                theme_health_disease_2,
                                width_two_plots #c()
                                
){
  list_plot <- list()
  tobind <- list()
  health_disease <- list()
  list_legend <- list()
  
  for (dataframes in names(final_df)){
    
    data <- final_df[[dataframes]]
    
    if (dataframes %in% names_df_tobind){
      tobind[[dataframes]] <- format_dataframe[[dataframes]]
      colnames(tobind[[dataframes]]) <- colnames(data)
      tobind[[dataframes]] <- tobind_row(tobind[[dataframes]], rows_tobind[[dataframes]])
      tobind[[dataframes]] <- tobind[[dataframes]] %>%
        mutate_at(vars(all_of(numeric_columns)), as.numeric)
    }
    
    
    if (dataframes %in% names_df_tobind){
      
      data <- rbind(data, tobind[[dataframes]])
    }
    
    data <- data %>%
      filter(!(!!sym(clinical_var_col) %in% to_exclude)) %>%
      group_by(!!sym(clinical_var_col)) %>%
      mutate(sum_ = sum(!!sym(pos_val)),
             sum_n = sum(!!sym(neg_val)) * -1) %>%
      ungroup()
    
    
    data <- data %>%
      mutate(sum_pos = sum_ + asterisk_pos[[dataframes]]) %>%
      mutate(sum_neg = sum_n - asterisk_neg[[dataframes]]) %>%
      select(-c(sum_, sum_n))
    
    data <- data %>%
      mutate(sum_arrange = !!sym(pos_val) + !!sym(neg_val)) %>%
      ungroup()
    
    data[[neg_val]] <- data[[neg_val]] * -1
    
    databis <- data %>%
      pivot_longer(
        cols = starts_with("correlation"),
        names_to = name_type_cor_col,
        values_to = name_type_cor_val) %>%
      mutate(!!sym(name_type_cor_col) := case_when(
        !!sym(name_type_cor_col) == pos_val ~ rename_pos_val,
        !!sym(name_type_cor_col) == neg_val ~ rename_neg_val,
        TRUE ~ !!sym(name_type_cor_col)))
    
    databis$color <- paste0(databis$type, databis[[name_type_cor_col]])
    
    to_add_y <- ifelse_condition_y(dataframes)
    
    
    databis <- databis %>%
      group_by(!!sym(clinical_var_col), type) %>%
      mutate(chi2 = ifelse(!!sym(name_type_cor_val) == min(!!sym(name_type_cor_val)), 1, chi2))
    
    
    databis <- databis %>%
      mutate(asterisk = case_when(
        chi2 > 0.01 & chi2 < 0.05 ~ "*",
        chi2 > 0.001 & chi2 <= 0.01 ~ "**",
        chi2 > 0.0001 & chi2 <= 0.001 ~ "***",
        chi2 <= 0.0001 ~ "****",
        TRUE ~ ""
      ))
    
    
    databis <- databis %>%
      mutate(!!sym(clinical_var_col) := gsub("_", " ", as.character(!!sym(clinical_var_col))))
    
    
    databis <- databis %>%
      arrange(clinical_var_col, sum_arrange) %>%
      select(-sum_arrange)
    
    
    databis_enriched <- databis %>%
      filter(type == "enriched")
    
    databis_depleted <- databis %>%
      filter(type == "depleted")
    
    
    for (which_plot in which_kind_plot){
      
      if (which_plot == "health"){
        
        clinical_filter <- list_clinical_variable[[1]]
        length_clin <- length(list_clinical_variable[[1]])
        
      } else if (which_plot == "disease"){
        
        
        clinical_filter  <- list_clinical_variable[[2]]
        length_clin <- length(list_clinical_variable[[2]])
        
      } else {
        
        clinical_filter <- c(list_clinical_variable[[1]], list_clinical_variable[[2]])
        length_clin <- length(c(list_clinical_variable[[1]],list_clinical_variable[[2]]))
      }
      
      
      asterisk_txt_pos <- geom_text(data = databis_enriched %>%
                                      filter(!!sym(clinical_var_col) %in% clinical_filter),
                                    aes(x = !!sym(clinical_var_col),
                                        y = (sum_pos),
                                        label = asterisk),
                                    vjust = -0.5, size = size_asterisk[[dataframes]], color = "black", angle = 0)
      
      asterisk_txt_neg <- geom_text(data = databis_depleted %>%
                                      filter(!!sym(clinical_var_col) %in% clinical_filter),
                                    aes(x = !!sym(clinical_var_col),
                                        y = sum_neg,
                                        label = asterisk),
                                    vjust = -0.5, size = size_asterisk[[dataframes]], color = "black", angle = 0)
      
      
      health_disease[[length(health_disease) +1]] <-
        ggplot(databis %>%
                 filter(!!sym(clinical_var_col) %in% clinical_filter) %>%
                 mutate(!!sym(clinical_var_col) := factor(!!sym(clinical_var_col), level = clinical_filter)) %>%
                 mutate(color = factor(color, level = c("enrichedpositive", "enrichednegative",
                                                        "depletedpositive", "depletednegative"))),
               aes(x= !!sym(clinical_var_col))) +
        
        
        geom_bar(aes(y=ifelse(correlation_value > 0, correlation_value, 0), fill = color), stat = "identity", color = "black") +
        geom_bar(aes(y=ifelse(correlation_value < 0, correlation_value, 0), fill = color), stat = "identity", color = "black") +
        
        theme_classic() +
        
        
        scale_y_continuous(breaks = limits_by[[dataframes]] ,
                           limits = limits_[[dataframes]],
                           labels = limits_by[[dataframes]]) +
        
        labs(y = paste(number_of[1], to_add_y, number_of[2], sep = " ")) +
        
        labs_plot[[dataframes]] +
        
        
        theme_size[[dataframes]]  +
        
        
        scale_fill_manual(values = colors_fill,
                          labels = labels_legend,
                          breaks = keep_in_legend) +
        
        geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        
        
        asterisk_txt_pos +
        
        asterisk_txt_neg
      
    }
    
    if (dataframes == "HC_filter"){
      
      gb_bis_ <- ggplot_gtable(ggplot_build(health_disease[[2]]))
      leg_bis_ <- gb_bis_$grobs[which(sapply(gb_bis_$grobs, function(x) x$name) == "guide-box")]
      leg_bis_ <- leg_bis_[[1]]
      leg_bis_$layout[leg_bis_$layout$name == "guide-box"]$height <- unit(0, "pt")
      list_legend[[length(list_legend) +1]] <- leg_bis_
      
    }
    
    
    
    p <- ggarrange(
      health_disease[[1]] + theme_health_disease_1[[dataframes]]  + coord_cartesian(clip = "off"),
      health_disease[[2]] + theme_health_disease_2[[dataframes]]  + coord_cartesian(clip = "off"),
      ncol = 2,
      widths = width_two_plots
    )
    
    health_disease <- list()
    
    p2 <- p +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      
      annotate("text", x = x_legend[[dataframes]][1], y =y_legend[[dataframes]][1], label =  label_correlation[1], size = size_legend[[dataframes]][1], fontface = "italic") +
      
      annotate("text", x = x_legend[[dataframes]][2], y = y_legend[[dataframes]][2], label =  label_correlation[2], size = size_legend[[dataframes]][1], fontface = "italic")
    
    
    p2 <- p2 +
      labs(title = letter_fig[[dataframes]]) +
      theme_plot_title[[dataframes]]
    
    
    list_plot[[length(list_plot) +1]] <- p2
    
  }
  
  
  return(list(plots = list_plot, legends = list_legend))
}

library(circlize)
## Function use for Heatmap (Figure)

heatmap_plot <- function(mat_cor, 
                         mat_corP, 
                         lgd_title = NULL,
                         size_lgd_gp = NULL ,
                         size_label_gp = NULL,
                         rownames_size = NULL,
                         colnames_size = NULL,
                         cluster_rows = TRUE, 
                         annot = NULL, 
                         cutoff_filled = 0.001) {
  
  col_fun <- colorRamp2(c(min(mat_cor), 0, max(mat_cor)), c("#D60C00FF", "white", "#072AC8")) 
  
  ht <- Heatmap(
    as.matrix(mat_cor), 
    cluster_rows = cluster_rows, 
    cluster_columns = TRUE, 
    col = col_fun,
    right_annotation = annot,
    row_names_gp = gpar(fontsize = rownames_size, col = "black"),
    column_names_gp = gpar(fontsize = colnames_size, col = "black"),
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(mat_corP[i, j])) {
        if (mat_corP[i, j] < cutoff_filled) {
          grid.points(pch = 19, x, y, size = unit(0.5, "char"))
        }
        if (mat_corP[i, j] < 0.05 & mat_corP[i, j] >= cutoff_filled) {
          grid.points(pch = 1, x, y, size = unit(0.5, "char")) 
        }
      }
      
      grid.rect(x = x, y = y, width = w, height = h, 
                gp = gpar(col = "black", fill = NA))  
    }
  )
  return(ht)
}



