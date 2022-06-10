#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Predicts HRS/WLS/UKB phenotypes from scores created using public data
# Date: 07/13/2020
# Author: Joel Becker & Aysu Okbay 

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "stringr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# Parse arguments
args=commandArgs(trailingOnly=TRUE)
cohort=args[1]

# Score types: Public, single-trait, multi-trait
#score_types=c("public","single","multi")
score_types=c("single_SBayesR","single_LDpred")


########################################################
###################### Data paths ######################
########################################################

# Phenotype directory
pheno_wd <- paste0("derived_data/10_Prediction/input/",cohort)

# Phenotype files
pheno_files <- list.files(pheno_wd)

# Phenotype names
pheno_names <- gsub(".pheno", "", pheno_files)

predict <- c()

# Loop over PGI types
for ( i in score_types ){
    # PGI directory
    assign(paste0(i,"_score_wd"), paste0("derived_data/9_Scores/",i,"/scores/"))
    assign(paste0(i,"_score_files"), list.files(eval(parse(text=paste0(i,"_score_wd"))), cohort))
    
    # PGI names
    assign(paste0(i,"_score_names"), gsub(paste0("PGS_",cohort,"_"), "", gsub("-.*", "",eval(parse(text=paste0(i,"_score_files"))))))
  
    # Overlapping names (both pheno and PGI available)
    assign(paste0("predict_",i), pheno_names[pheno_names %in% eval(parse(text=paste0(i,"_score_names")))])
    
    predict <- union(predict, eval(parse(text=paste0("predict_",i))))
}
pheno_names <- predict


# If cohort is HRS, need pheno-geno crosswalk
if ( cohort == "HRS2" ){
    # Scores-phenos crosswalk
    score_pheno_crosswalk_path <- paste0("original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.dta")
    score_pheno_crosswalk_data <- read.dta(score_pheno_crosswalk_path) %>%
        mutate(IID = as.numeric(LOCAL_ID),
        HHID = as.numeric(HHID),
        PN = as.numeric(PN))
}

# PCs
if ( cohort == "UKB3" ){
    PCs_path <- "derived_data/10_Prediction/input/UKB3/PC_BATCHdum.txt"
    PCs_data <- fread(PCs_path)
} else {
    PCs_path <- paste0("derived_data/8_PCs/",cohort,"/",cohort,"_PCs.eigenvec")
    PCs_oldnames <- paste0("V", 3:22)
    PCs_newnames <- paste0("pc", 1:20)
    PCs_data <- fread(PCs_path) %>%
    rename(IID = V2) %>%
    rename_at(vars(PCs_oldnames), ~ PCs_newnames)
}

### Empty data frame with r^2
df <- data.frame(phenotype = character())
bootstrap_df <- data.frame(
  iterations = numeric(),
  phenotype = character(),
  inc_r2_public_LDpred = numeric(),
  inc_r2_single_LDpred = numeric(),
  inc_r2_multi_LDpred = numeric(),
  inc_r2_single_SBayesR = numeric()
)
#####################################################3

# Set number of iterations for bootstrapping
iterations <- 1000

# Loop over phenotypes
for (i in 1:length(pheno_names)){
    
    local_df <- data.frame(phenotype = character(1))
    pheno <- pheno_names[i]
    local_df$phenotype <- pheno

    # Read in pheno data
    pheno_path <- paste0(pheno_wd,"/", pheno, ".pheno")
    pheno_data <- fread(pheno_path)
    
    score_data <- NULL
    # Loop over score types
    for (i in score_types) {
        score_files <- eval(parse(text=paste0(i,"_score_files")))   
        
        # Score file for pheno
        score_file <- score_files[grep(paste0("_",pheno,"-"), score_files, fixed=TRUE)]

        if ( length(score_file) != 0 ) {
            # Dummy for whether score exists
            assign(paste0(pheno,"_",i),"TRUE")

            # Score path
            score_path <- paste0(eval(parse(text=paste0(i,"_score_wd"))),score_file)

            # Read in score data
            local_score_data <- fread(score_path)
            names(local_score_data)[5] <- paste0("score_",i)
    
            # Merge with other score types
            if ( is.null(score_data) ) {
                score_data <- local_score_data
            } else {
                score_data <- inner_join(score_data, local_score_data, by="IID")
            }
        } else {
            assign(paste0(pheno,"_",i),"FALSE")
        }
    }  
    
    # Merge with pheno and PC data
    # HRS
    if ( cohort == "HRS2" ){
        if ("hhidpn" %in% colnames(pheno_data)){ # if needing to join on hhidpn variables
            score_pheno_crosswalk_data <- score_pheno_crosswalk_data %>% mutate(hhidpn = (1000 * HHID) + PN)
        
            data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
                inner_join(pheno_data, by="hhidpn") %>%
                inner_join(PCs_data, by="IID") %>%
                drop_na(phenotype,pc1)
            
        } else { # else, if joining on hhid and pn
            if ("hhid" %in% colnames(pheno_data)) {
                pheno_data <- pheno_data %>% mutate(HHID = hhid, PN = pn)
            }
            data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
                inner_join(pheno_data, by=c("HHID", "PN")) %>%
                inner_join(PCs_data, by="IID") %>%
                drop_na(phenotype,pc1)
        }
    } else if ( cohort == "WLS" ){
        # WLS
        score_data <- score_data %>%
            mutate(id = as.numeric(IID))
        data <- pheno_data %>%
            inner_join(score_data, by="id") %>%
            inner_join(PCs_data, by="IID") %>%
            drop_na(phenotype,pc1)
    } else {
        # UKB3
        data <- pheno_data %>%
            inner_join(score_data, by="IID") %>%
            inner_join(PCs_data, by="IID") %>%
            drop_na(phenotype,PC1)
    }

    # Regress with and without scores, extract r^2
    PCs <- str_c(paste0("pc", 1:10), collapse=" + ")
    if ( cohort == "UKB3" ){
        PCs <- str_c(paste0("PC",1:20),collapse=" + ")
        BATCH <- str_c(paste0("batch",1:106),collapse=" + ")
        PCs <- paste(PCs,BATCH,sep=" + ")
    }
    
    formula_no_score <- as.formula(paste0("phenotype ~ ", PCs))
    reg_no_score <- lm(formula_no_score, data)
    local_df$r2_no_score <- summary(reg_no_score)$r.squared * 100

    local_df$r2_inc_single_LDpred <- NA
    local_df$r2_inc_single_SBayesR <- NA
    #local_df$r2_inc_multi_LDpred <- NA
    #local_df$r2_inc_public_LDpred <- NA
    
    
    for (k in score_types){
        assign(paste0("formula_with_",k,"_score"),as.formula(paste0("phenotype ~ score_", k, " + ", PCs)))
        reg <- eval(parse(text=paste0(pheno,"_",k)))
            if ( reg == TRUE ){
                reg_with_score <- lm(eval(parse(text=paste0("formula_with_",k,"_score"))), data)
                r2 <-summary(reg_with_score)$r.squared*100
            } else {
                r2 <- NA
            }
        local_df[1,paste0("r2_inc_",k)] <- r2 - local_df$r2_no_score[1]
    }

    #local_df$r2_diff_single_LDpred_public_LDpred  <- local_df$r2_inc_single_LDpred - local_df$r2_inc_public_LDpred
    #local_df$r2_diff_multi_LDpred_public_LDpred  <- local_df$r2_inc_multi_LDpred - local_df$r2_inc_public_LDpred
    #local_df$r2_diff_multi_LDpred_single_LDpred  <- local_df$r2_inc_multi_LDpred - local_df$r2_inc_single_LDpred
    local_df$r2_diff_single_SBayesR_single_LDpred  <- local_df$r2_inc_single_SBayesR - local_df$r2_inc_single_LDpred
    
    # Vectors for bootstrapping data
    for (k in score_types){
        assign(paste0("r2_inc_",k,"_list"), c())
    }
    #r2_diff_single_LDpred_public_LDpred_list <- c()
    #r2_diff_multi_LDpred_public_LDpred_list <- c()
    #r2_diff_multi_LDpred_single_LDpred_list <- c()
    r2_diff_single_SBayesR_single_LDpred_list <-c()

    # Start bootstrapping R2's
    for (j in 1:iterations) {
        set.seed(j)
        
        # print iteration and phenotype
        if (j %% 100 == 0){
            print(paste0("Iteration: ", j, ". Phenotype: ", pheno, "."))
        }

        # resample, regress, add r2
        resampled_data <- sample_n(data, nrow(data), replace=TRUE)
        resampled_no_score_reg <- lm(formula_no_score, resampled_data)
        resampled_r2_no_score <- summary(resampled_no_score_reg)$r.squared * 100
        
        for (k in score_types){
            reg=eval(parse(text=paste0(pheno,"_",k)))
            if ( reg == TRUE ){
                resampled_reg_with_score <- lm(eval(parse(text=paste0("formula_with_",k,"_score"))), resampled_data)
                resampled_r2 <-summary(resampled_reg_with_score)$r.squared*100
            } else {
                resampled_r2 <- NA
            }

            assign(eval(paste0("r2_inc_",k,"_list")), c(eval(parse(text=paste0("r2_inc_",k,"_list"))), resampled_r2-resampled_r2_no_score))
        }    
        
    }
    r2_diff_single_SBayesR_single_LDpred_list <- r2_inc_single_SBayesR_list - r2_inc_single_LDpred_list
    #r2_diff_single_LDpred_public_LDpred_list <- r2_inc_single_LDpred_list - r2_inc_public_LDpred_list
    #r2_diff_multi_LDpred_public_LDpred_list <- r2_inc_multi_LDpred_list - r2_inc_public_LDpred_list
    #r2_diff_multi_LDpred_single_LDpred_list <- r2_inc_multi_LDpred_list - r2_inc_single_LDpred_list 

    # save bootstrapped r2s
    local_bootstrap_df <- data.frame(
        iteration = 1:iterations,
        phenotype = pheno
    )
    bootstrap_df <- rbind(bootstrap_df, local_bootstrap_df)

    for (k in score_types){
        assign(eval(paste0("local_bootstrap_df$inc_r2_",k)), eval(parse(text=paste0("r2_inc_",k,"_list"))))
        assign(paste0("r2_inc_",k,"_low"), quantile(eval(parse(text=paste0("r2_inc_",k,"_list"))),.025, na.rm=T)[[1]])
        assign(paste0("r2_inc_",k,"_high"), quantile(eval(parse(text=paste0("r2_inc_",k,"_list"))),.975, na.rm=T)[[1]])
        assign(paste0("r2_inc_",k,"_mean"), mean(eval(parse(text=paste0("r2_inc_",k,"_list")))))
    }

    r2_diff_single_SBayesR_single_LDpred_low  <- quantile(r2_diff_single_SBayesR_single_LDpred_list, .025, na.rm=T)[[1]]
    r2_diff_single_SBayesR_single_LDpred_high  <- quantile(r2_diff_single_SBayesR_single_LDpred_list, .975, na.rm=T)[[1]]
    r2_diff_single_SBayesR_single_LDpred_mean <- mean(r2_diff_single_SBayesR_single_LDpred_list)
    
    #r2_diff_single_LDpred_public_LDpred_low  <- quantile(r2_diff_single_LDpred_public_LDpred_list, .025, na.rm=T)[[1]]
    #r2_diff_single_LDpred_public_LDpred_high  <- quantile(r2_diff_single_LDpred_public_LDpred_list, .975, na.rm=T)[[1]]
    #r2_diff_single_LDpred_public_LDpred_mean <- mean(r2_diff_single_LDpred_public_LDpred_list)

    #r2_diff_multi_LDpred_public_LDpred_low  <- quantile(r2_diff_multi_LDpred_public_LDpred_list, .025, na.rm=T)[[1]]
    #r2_diff_multi_LDpred_public_LDpred_high  <- quantile(r2_diff_multi_LDpred_public_LDpred_list, .975, na.rm=T)[[1]]
    #r2_diff_multi_LDpred_public_LDpred_mean <- mean(r2_diff_multi_LDpred_public_LDpred_list)

    #r2_diff_multi_LDpred_single_LDpred_low  <- quantile(r2_diff_multi_LDpred_single_LDpred_list, .025, na.rm=T)[[1]]
    #r2_diff_multi_LDpred_single_LDpred_high  <- quantile(r2_diff_multi_LDpred_single_LDpred_list, .975, na.rm=T)[[1]]
    #r2_diff_multi_LDpred_single_LDpred_mean <- mean(r2_diff_multi_LDpred_single_LDpred_list)

    
    # print r2 confidence interval
    print(paste0("Incremental R^2 single LDpred(95% CI): ", r2_inc_single_LDpred_mean, " [", r2_inc_single_LDpred_low, " : ", r2_inc_single_LDpred_high, "]"))
    print(paste0("Incremental R^2 single SBayesR(95% CI): ", r2_inc_single_SBayesR_mean, " [", r2_inc_single_SBayesR_low, " : ", r2_inc_single_SBayesR_high, "]"))
    print(paste0("Diff R^2 single_SBayesR - single LDpred(95% CI): ", r2_diff_single_SBayesR_single_LDpred_mean, " [", r2_diff_single_SBayesR_single_LDpred_low, " : ", r2_diff_single_SBayesR_single_LDpred_high, "]"))
    
    #print(paste0("Incremental R^2 public LDpred(95% CI): ", r2_inc_public_LDpred_mean, " [", r2_inc_public_LDpred_low, " : ", r2_inc_public_LDpred_high, "]"))
    #print(paste0("Incremental R^2 multi LDpred(95% CI): ", r2_inc_multi_LDpred_mean, " [", r2_inc_multi_LDpred_low, " : ", r2_inc_multi_LDpred_high, "]"))
    #print(paste0("Diff R^2 single LDpred - public LDpred(95% CI): ", r2_single_LDpred_public_LDpred_mean, " [", r2_single_LDpred_public_LDpred_low, " : ", r2_diff_single_public_high, "]"))
    #print(paste0("Diff R^2 multi LDpred - public LDpred(95% CI): ", r2_diff_multi_LDpred_public_LDpred_mean, " [", r2_diff_multi_LDpred_public_LDpred_low, " : ", r2_diff_multi_LDpred_public_LDpred_high, "]"))
    #print(paste0("Diff R^2 multi LDpred - single LDpred(95% CI): ", r2_diff_multi_LDpred_single_LDpred_mean, " [", r2_diff_multi_LDpred_single_LDpred_low, " : ", r2_diff_multi_LDpred_single_LDpred_high, "]"))
    

    for (k in score_types){
        assign(eval(paste0("local_df$r2_inc_",k,"_lower")), eval(parse(text=paste0("r2_inc_",k,"_low"))))
        assign(eval(paste0("local_df$r2_inc_",k,"_upper")), eval(parse(text=paste0("r2_inc_",k,"_high"))))
    }
    #local_df$r2_diff_single_LDpred_public_LDpred_lower <- r2_diff_single_LDpred_public_LDpred_low
    #local_df$r2_diff_single_LDpred_public_LDpred_upper <- r2_diff_single_LDpred_public_LDpred_high

    #local_df$r2_diff_multi_LDpred_public_LDpred_lower <- r2_diff_multi_LDpred_public_LDpred_low
    #local_df$r2_diff_multi_LDpred_public_LDpred_upper <- r2_diff_multi_LDpred_public_LDpred_high

    #local_df$r2_diff_multi_LDpred_single_LDpred_lower <- r2_diff_multi_LDpred_single_LDpred_low
    #local_df$r2_diff_multi_LDpred_single_LDpred_upper <- r2_diff_multi_LDpred_single_LDpred_high

    local_df$r2_diff_single_SBayesR_single_LDpred_lower <- r2_diff_single_SBayesR_single_LDpred_low
    local_df$r2_diff_single_SBayesR_single_LDpred_upper <- r2_diff_single_SBayesR_single_LDpred_high
    
    # print and record merged prediction sample N
    local_df$N <- nrow(data)
    print(paste0("N: ", nrow(data)))

    # bind local df to master df
    df <- rbind(local_df, df)
}


########################################################
#################### Save results ######################
########################################################

# low-level bootstraps
fwrite(
  bootstrap_df,
  paste0("derived_data/10_Prediction/output/bootstraps/",cohort,"_bootstraps.txt")
)

#results
fwrite(
  df,
  paste0("derived_data/10_Prediction/output/",cohort,"_r2.txt")
)
