#!/usr/bin/env Rscript


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "stringr","ggplot2","gridExtra")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# Parse arguments
args=commandArgs(trailingOnly=TRUE)
cohort=args[1]
score_type=args[2]
mainDir=args[3]


########################################################

# PCs
if ( cohort == "UKB3" | cohort == "UKB2" | cohort == "UKB1" ){
    PCs_path <- paste0(mainDir,"/derived_data/10_Prediction/input/UKB3/PC_BATCHdum.txt")
    PCs_data <- fread(PCs_path)
} else {
    PCs_path <- paste0(mainDir,"/derived_data/8_PCs/",cohort,"/",cohort,"_PCs.eigenvec")
    PCs_oldnames <- paste0("V", 3:22)
    PCs_newnames <- paste0("pc", 1:20)
    PCs_data <- fread(PCs_path) %>%
    rename(IID = V2) %>%
    rename_at(vars(PCs_oldnames), ~ PCs_newnames)
}

score_wd <- paste0(mainDir,"/derived_data/9_Scores/",score_type,"/scores/")

score_files <- list.files(score_wd)
score_files <- score_files[grep(cohort, score_files)]
score_names <- gsub(paste0("PGS_",cohort,"_"), "", 
                    gsub("_LDpred_p1.txt", "", 
                            gsub("-.*", "",score_files)))

# Merge scores
score_data <- NULL
for (i in 1:length(score_names)){
    score <- score_names[i]
    score_file <- score_files[grep(paste0("_",score,"-"), score_files, fixed=TRUE)]
    score_path <- paste0(score_wd,score_file)

    # Read in score data
    local_score_data <- fread(score_path)
    names(local_score_data)[5] <- score_names[i]
    
    if ( is.null(score_data) ) {
        score_data <- local_score_data
    } else {
        score_data <- inner_join(score_data, local_score_data, by=c("FID","IID"))
    }
}        
score_data <- score_data %>% select(-contains("ALLELE"))

# Merge with PC data
if ( cohort == "UKB3" | cohort == "UKB2" | cohort == "UKB1" ){
    data <- score_data %>%
        inner_join(PCs_data, by="IID") #%>%
        drop_na(PC1,EA)
        print("4")
} else {
    data <- score_data %>%
        inner_join(PCs_data, by="IID") %>%
        drop_na(pc1,EA)
}

# Residualize on PCs
PCs <- str_c(paste0("pc", 1:10), collapse=" + ")
if ( cohort == "UKB3" | cohort == "UKB2" | cohort == "UKB1" ){
        PCs <- str_c(paste0("PC",1:20),collapse=" + ")
        BATCH <- str_c(paste0("batch",1:106),collapse=" + ")
        PCs <- paste(PCs,BATCH,sep=" + ")
}

pdf(file=paste0("Histograms_",cohort,"_",score_type,".pdf"))
for (i in 1:length(score_names)){
    score <- score_names[i]
    print(score)
    # Histogram w/o residualization        
    hist <- qplot(data[[score]], 
        geom = "histogram",  
        main = paste0(score,"- not residualized"), 
        xlab = paste0("PGI - ",score),
        fill=I("grey"), 
        col=I("black"),
        bins=50)
    model <- as.formula(paste0(score," ~ ", PCs))
    reg <- summary(lm(model, data))
    hist_resid <- qplot(reg$resid, 
        geom="histogram",  
        main = paste0(score,"- residualized"), 
        xlab = paste0("PGI - ",score),
        fill=I("grey"), 
        col=I("black"),
        bins=50)
    grid.arrange(hist, hist_resid, nrow = 2, newpage = F)
    grid::grid.newpage()
}

dev.off()
