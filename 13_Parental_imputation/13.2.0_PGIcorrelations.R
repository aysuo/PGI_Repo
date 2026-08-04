
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
pgi_dir=args[2]
out_dir=args[3]

########################################################
###################### Data paths ######################
########################################################

parental_PGI_files <- list.files(pgi_dir, paste0(cohort, "_parental"))
phenos <- gsub(paste0("PGS_",cohort,"_parental_"), "", gsub("-.*", "",parental_PGI_files))

corResults <- data.frame(phenotype = character(length(phenos)))

# Loop over phenotypes
for (i in 1:length(phenos)){
    
    pheno <- phenos[i]
    corResults$phenotype[i] <- pheno

    # Read in parental PGI file
    parentalPGI <- fread(paste0(pgi_dir,"/",parental_PGI_files[i]))
    regularPGI <- fread(paste0(pgi_dir,"/PGS_", cohort , "_" , pheno ,"-single_SBayesR.txt"))

    PGI <- inner_join(parentalPGI, regularPGI, by="IID")

    if ( "maternal" %in% colnames(PGI) ){
        corResults$proband_maternal[i]  <- cor(PGI$proband, PGI$maternal)
        corResults$proband_paternal[i]  <- cor(PGI$proband, PGI$paternal)
        corResults$proband_parental[i]  <- NA
    }
    else {
        corResults$proband_maternal[i]  <- NA
        corResults$proband_paternal[i]  <- NA
        corResults$proband_parental[i]  <- cor(PGI$proband, PGI$parental)
    } 

    corResults$probandRegular_probandSnipar[i] <- cor(PGI$proband, PGI[,paste0("PGS_", cohort, "_", pheno,"-single_SBayesR"),with=FALSE])

    corResults$N[i]  <- nrow(PGI)
}
    
fwrite(
  corResults,
  paste0(out_dir,"/",cohort,"_cors.txt"),
  sep="\t",
  quote=F
)

