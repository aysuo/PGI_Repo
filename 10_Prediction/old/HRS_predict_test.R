#----------------------------------------------------------------------------------#
# Predicts HRS phenotypes from scores created using public data
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

# Set directory
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction")


########################################################
###################### Data paths ######################
########################################################

# phenos
pheno_wd <- "input/HRS"
pheno_files <- list.files(pheno_wd)
pheno_names <- gsub(".pheno", "", pheno_files)

score_types=c("public","single","multi")

for(i in score_types){
    assign(paste0(i,"_score_wd"), paste0("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/",i,"/scores/"))
    assign(paste0(i,"_score_files"), list.files(eval(parse(text=paste0(i,"_score_wd")))))
    
    # score files
    score_files <- eval(parse(text=paste0(i,"_score_files")))

    # HRS score files
    assign(paste0("HRS_",i,"_score_files"),score_files[grep("HRS", score_files)])

    # HRS score names
    assign(paste0(i,"_score_names"), gsub("PGS_HRS2_", "", gsub("_LDpred_p1.txt", "", 
        gsub("-.*", "",eval(parse(text=paste0("HRS_",i,"_score_files")))))))

    # Overlapping names (both pheno and score)
    assign(paste0("predict_",i), pheno_names[pheno_names %in% eval(parse(text=paste0(i,"_score_names")))])
}

# List of phenotypes for which there's at least one score
pheno_names <- union(predict_single, predict_multi)

# Scores-phenos crosswalk
score_pheno_crosswalk_path <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.dta"
score_pheno_crosswalk_data <- read.dta(score_pheno_crosswalk_path) %>%
  mutate(IID = as.numeric(LOCAL_ID),
         HHID = as.numeric(HHID),
         PN = as.numeric(PN))

# PCs
PCs_path <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/8_PCs/HRS2/HRS2_PCs.eigenvec"
PCs_oldnames <- paste0("V", 3:22)
PCs_newnames <- paste0("pc", 1:20)
PCs_data <- fread(PCs_path) %>%
  rename(IID = V2) %>%
  rename_at(vars(PCs_oldnames), ~ PCs_newnames)

### empty data frame with r^2
df <- data.frame(phenotype = character())
bootstrap_df <- data.frame(
  iterations = numeric(),
  phenotype = character(),
  inc_r2_public = numeric(),
  inc_r2_single = numeric(),
  inc_r2_multi = numeric()
)
#####################################################3
# options for prediction loop
iterations <- 1000

for (i in 1:length(pheno_names)){

    local_df <- data.frame(phenotype = character(1))
    pheno <- pheno_names[i]
    local_df$phenotype <- pheno

    # Read in pheno data
    pheno_path <- paste0(getwd(), "/input/HRS/", pheno, ".pheno")
    pheno_data <- fread(pheno_path)

    score_data <- NULL
    for (i in score_types) {
        score_files <- eval(parse(text=paste0("HRS_",i,"_score_files")))   

        # Score file for pheno
        score_file <- score_files[grep(paste0("_",pheno), score_files, fixed=TRUE)]

        if ( length(score_file) != 0 ) {
            # Dummy for whether score exists
            assign(paste0(pheno,"_",i),"TRUE")

            # Score path
            score_path <- paste0(eval(parse(text=paste0(i,"_score_wd"))),score_file)

            # Read in score data
            local_score_data <- fread(score_path)
            names(local_score_data)[5] <- paste0("score_",i)
    
            if ( is.null(score_data) ) {
                score_data <- local_score_data
            } else {
                score_data <- inner_join(score_data, local_score_data, by="IID")
            }
        } else {
            assign(paste0(pheno,"_",i),"FALSE")
        }
    }  

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

    # regress with and without scores, extract r^2
    PCs <- str_c(paste0("pc", 1:10), collapse=" + ")
    formula_no_score <- as.formula(paste0("phenotype ~ ", PCs))
    formula_with_public_score <- as.formula(paste0("phenotype ~ score_public + ", PCs))
    formula_with_single_score <- as.formula(paste0("phenotype ~ score_single + ", PCs))
    formula_with_multi_score <- as.formula(paste0("phenotype ~ score_multi + ", PCs))
    
    reg_no_score <- lm(formula_no_score, data)
    local_df$r2_no_score <- summary(reg_no_score)$r.squared * 100

    for (k in score_types){
        reg <- eval(parse(text=paste0(pheno,"_",k)))
            if ( reg == TRUE ){
                assign(paste0("reg_with_",k,"_score"), lm(eval(parse(text=paste0("formula_with_",k,"_score"))), data))
                assign(paste0("r2_with_",k,"_score"), eval(parse(text=paste0("summary(reg_with_",k,"_score)$r.squared"))) * 100)
            } else {
                assign(paste0("r2_with_",k,"_score"),NA)
            }
    }

    local_df$r2_with_public_score <- r2_with_public_score
    local_df$r2_with_single_score <- r2_with_single_score
    local_df$r2_with_multi_score <- r2_with_multi_score

    local_df$r2_inc_public <- local_df$r2_with_public_score - local_df$r2_no_score
    local_df$r2_inc_single <- local_df$r2_with_single_score - local_df$r2_no_score
    local_df$r2_inc_multi  <- local_df$r2_with_multi_score - local_df$r2_no_score

    local_df$r2_diff_single_public  <- local_df$r2_with_single_score - local_df$r2_with_public_score
    local_df$r2_diff_multi_public  <- local_df$r2_with_multi_score - local_df$r2_with_public_score
    local_df$r2_diff_multi_single  <- local_df$r2_with_multi_score - local_df$r2_with_single_score

    # bootstrap
    r2_inc_public_list <- c()
    r2_inc_single_list <- c()
    r2_inc_multi_list <- c()
    r2_diff_single_public_list <- c()
    r2_diff_multi_public_list <- c()
    r2_diff_multi_single_list <- c()


    for (j in 1:iterations) {
        # print iteration and phenotype
        if (j %% 100 == 0){
            print(paste0("Iteration: ", j, ". Phenotype: ", pheno, "."))
        }

        # resample, regress, add r2
        resampled_data <- sample_n(data, nrow(data), replace=TRUE)
        resampled_no_score_reg_fit <- lm(formula_no_score, resampled_data)
        resampled_no_score_reg <- summary(resampled_no_score_reg_fit)
        r2_no_score <- resampled_no_score_reg$r.squared * 100
        
        for (k in score_types){
            reg=eval(parse(text=paste0(pheno,"_",k)))
            if ( reg == TRUE ){
                assign(paste0("resampled_with_",k,"_score_reg_fit"), lm(eval(parse(text=paste0("formula_with_",k,"_score"))), resampled_data))
                assign(paste0("r2_with_",k,"_score"), eval(parse(text=paste0("summary(resampled_with_",k,"_score_reg_fit)$r.squared"))) * 100)
            } else {
                assign(paste0("r2_with_",k,"_score"), NA)
            }
        }

        r2_inc_public_list[j] <- r2_with_public_score - r2_no_score
        r2_inc_single_list[j] <- r2_with_single_score - r2_no_score
        r2_inc_multi_list[j] <- r2_with_multi_score - r2_no_score

        r2_diff_single_public_list[j] <- r2_with_single_score - r2_with_public_score
        r2_diff_multi_public_list[j] <- r2_with_multi_score - r2_with_public_score
        r2_diff_multi_single_list[j] <- r2_with_multi_score - r2_with_single_score
        
    }

    # save bootstrapped r2s
    local_bootstrap_df <- data.frame(
        iteration = 1:1000,
        phenotype = pheno
    )
    

    local_bootstrap_df$inc_r2_public <- r2_inc_public_list
    local_bootstrap_df$inc_r2_single <- r2_inc_single_list
    local_bootstrap_df$inc_r2_multi <- r2_inc_multi_list

    bootstrap_df <- rbind(bootstrap_df, local_bootstrap_df)

    # generate quantities of interest
    r2_inc_public_low  <- quantile(r2_inc_public_list, .025, na.rm=T)[[1]]
    r2_inc_public_high <- quantile(r2_inc_public_list, .975, na.rm=T)[[1]]
    r2_inc_public_mean <- mean(r2_inc_public_list)

    r2_inc_single_low  <- quantile(r2_inc_single_list, .025, na.rm=T)[[1]]
    r2_inc_single_high <- quantile(r2_inc_single_list, .975, na.rm=T)[[1]]
    r2_inc_single_mean <- mean(r2_inc_single_list)

    r2_inc_multi_low  <- quantile(r2_inc_multi_list, .025, na.rm=T)[[1]]
    r2_inc_multi_high <- quantile(r2_inc_multi_list, .975, na.rm=T)[[1]]
    r2_inc_multi_mean <- mean(r2_inc_multi_list)

    r2_diff_single_public_low  <- quantile(r2_diff_single_public_list, .025, na.rm=T)[[1]]
    r2_diff_single_public_high  <- quantile(r2_diff_single_public_list, .975, na.rm=T)[[1]]
    r2_diff_single_public_mean <- mean(r2_diff_single_public_list)

    r2_diff_multi_public_low  <- quantile(r2_diff_multi_public_list, .025, na.rm=T)[[1]]
    r2_diff_multi_public_high  <- quantile(r2_diff_multi_public_list, .975, na.rm=T)[[1]]
    r2_diff_multi_public_mean <- mean(r2_diff_multi_public_list)

    r2_diff_multi_single_low  <- quantile(r2_diff_multi_single_list, .025, na.rm=T)[[1]]
    r2_diff_multi_single_high  <- quantile(r2_diff_multi_single_list, .975, na.rm=T)[[1]]
    r2_diff_multi_single_mean <- mean(r2_diff_multi_single_list)


    # print r2 confidence interval
    print(paste0("Incremental R^2 public (95% CI): ", r2_inc_public_mean, " [", r2_inc_public_low, " : ", r2_inc_public_high, "]"))
    print(paste0("Incremental R^2 single (95% CI): ", r2_inc_single_mean, " [", r2_inc_single_low, " : ", r2_inc_single_high, "]"))
    print(paste0("Incremental R^2 multi (95% CI): ", r2_inc_multi_mean, " [", r2_inc_multi_low, " : ", r2_inc_multi_high, "]"))
    print(paste0("Diff R^2 single - public (95% CI): ", r2_diff_single_public_mean, " [", r2_diff_single_public_low, " : ", r2_diff_single_public_high, "]"))
    print(paste0("Diff R^2 multi - public (95% CI): ", r2_diff_multi_public_mean, " [", r2_diff_multi_public_low, " : ", r2_diff_multi_public_high, "]"))
    print(paste0("Diff R^2 multi - single (95% CI): ", r2_diff_multi_single_mean, " [", r2_diff_multi_single_low, " : ", r2_diff_multi_single_high, "]"))

    ### bootstrap function
    local_df$r2_inc_public_lower <- r2_inc_public_low
    local_df$r2_inc_public_upper <- r2_inc_public_high

    local_df$r2_inc_single_lower <- r2_inc_single_low
    local_df$r2_inc_single_upper <- r2_inc_single_high

    local_df$r2_inc_multi_lower <- r2_inc_multi_low
    local_df$r2_inc_multi_upper <- r2_inc_multi_high

    local_df$r2_diff_single_public_lower <- r2_diff_single_public_low
    local_df$r2_diff_single_public_upper <- r2_diff_single_public_high

    local_df$r2_diff_multi_public_lower <- r2_diff_multi_public_low
    local_df$r2_diff_multi_public_upper <- r2_diff_multi_public_high

    local_df$r2_diff_multi_single_lower <- r2_diff_multi_single_low
    local_df$r2_diff_multi_single_upper <- r2_diff_multi_single_high

    
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

# overwrite main file
fwrite(
  bootstrap_df,
  paste0(
    getwd(),
    "/output/HRS_bootstraps.txt"
  )
)

# overwrite main file
fwrite(
  df,
  paste0(
    getwd(),
    "/output/HRS_phenotypes_r2.txt"
  )
)
