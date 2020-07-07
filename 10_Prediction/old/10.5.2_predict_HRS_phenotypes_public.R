#----------------------------------------------------------------------------------#
# Predicts HRS phenotypes from scores created using public data
# Date: 04/16/2020
# Author: Joel Becker

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
################ Negalkerke R^2 function ###############
########################################################

# please note:
# the below code is copied from https://github.com/cran/rcompanion/blob/master/R/nagelkerke.r
# (could not get the package to load on NBER servers)

nagelkerke <- function(fit, null=NULL, restrictNobs=FALSE) {

   TOGGLE =   (class(fit)[1]=="lm"
             | class(fit)[1]=="gls"
             | class(fit)[1]=="lme"
             | class(fit)[1]=="glm"
             | class(fit)[1]=="negbin"
             | class(fit)[1]=="zeroinfl"
             | class(fit)[1]=="clm"
             | class(fit)[1]=="vglm"
             | class(fit)[1]=="betareg"
             | class(fit)[1]=="rq")
   BOGGLE =   (class(fit)[1]=="nls"
             | class(fit)[1]=="lmerMod"
             | class(fit)[1]=="glmerMod"
             | class(fit)[1]=="merModLmerTest"
             | class(fit)[1]=="lmerModLmerTest"
             | class(fit)[1]=="clmm")
   SMOGGLE =   (class(fit)[1]=="lmerMod"
              | class(fit)[1]=="glmerMod"
              | class(fit)[1]=="merModLmerTest"
              | class(fit)[1]=="lmerModLmerTest"
              | class(fit)[1]=="vglm")
   ZOGGLE  = (class(fit)[1]=="zeroinfl")
   ZOGGLE2 = (class(fit)[1]=="rq")
   NOGGLE = is.null(null)
   ERROR  = "Note: For models fit with REML, these statistics are based on refitting with ML"
   ERROR2 = "None"

  if(!restrictNobs & NOGGLE  & TOGGLE){null = update(fit, ~ 1)}
  if(restrictNobs  & NOGGLE  & TOGGLE){null = update(fit, ~ 1, data=fit$model)}

  if(restrictNobs  & !NOGGLE){null = update(null, data=fit$model)}

  if(NOGGLE & BOGGLE)
     {ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"}
  if((!TOGGLE) & (!BOGGLE))
   {ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"}

   SMOGGLE2 = (class(null)[1]=="lmerMod"
              | class(null)[1]=="glmerMod"
              | class(null)[1]=="merModLmerTest"
              | class(null)[1]=="lmerModLmerTest"
              | class(null)[1]=="vglm")

  Y = matrix(rep(NA,2),
            ncol=1)
  colnames(Y) = ""
  rownames(Y) = c("Model:", "Null:")

  Z = matrix(rep(NA, 3),
             ncol=1)
  colnames(Z) = c("Pseudo.R.squared")
  rownames(Z) = c("McFadden", "Cox and Snell (ML)",
                  "Nagelkerke (Cragg and Uhler)")

  X = matrix(rep(NA,4),
             ncol=4)
  colnames(X) = c("Df.diff","LogLik.diff","Chisq","p.value")
  rownames(X) = ""

  U = matrix(rep(NA,2),
            ncol=1)
  colnames(U) = ""
  rownames(U) = c("Model:", "Null:")

  if(TOGGLE | BOGGLE){
  if (!SMOGGLE){Y[1]= toString(fit$call)}
  if (SMOGGLE){Y[1]= toString(fit@call)}
  }

  if(TOGGLE | (BOGGLE & !NOGGLE)){

  if (!SMOGGLE2){Y[2]= toString(null$call)}
  if (SMOGGLE2){Y[2]= toString(null@call)}

  if(!ZOGGLE & !ZOGGLE2){N = nobs(fit)
                         U[1,1]= nobs(fit); U[2,1]= nobs(null)}
  if(!ZOGGLE &  ZOGGLE2){N = length(fit$y)
                         U[1,1]= length(fit$y); U[2,1]= length(null$y)}
  if(ZOGGLE){N = fit$n
             U[1,1]= fit$n; U[2,1]= null$n}

  if(U[1,1] != U[2,1]){
    ERROR2 = "WARNING: Fitted and null models have different numbers of observations"}

  m = suppressWarnings(logLik(fit, REML=FALSE))[1]
  n = suppressWarnings(logLik(null, REML=FALSE))[1]
  mf = 1 - m/n
  Z[1,] = signif(mf, digits=6)
  cs = 1 - exp(-2/N * (m - n))
  Z[2,] = signif(cs, digits=6)
  nk = cs/(1 - exp(2/N * n))
  Z[3,] = signif(nk, digits=6)

  o = n - m
  dfm = attr(logLik(fit),"df")
  dfn = attr(logLik(null),"df")
  if(class(fit)[1]=="vglm"){dfm=df.residual(fit)}
  if(class(fit)[1]=="vglm"){dfn=df.residual(null)}
  dff = dfn - dfm
  CHI = 2 * (m - n)
  P = pchisq(CHI, abs(dff), lower.tail = FALSE)

  X [1,1] = dff
  X [1,2] = signif(o, digits=5)
  X [1,3] = signif(CHI, digits=5)
  X [1,4] = signif(P, digits=5)
  }

  W=ERROR

  WW=ERROR2

  V = list(Y, Z, X, U, W, WW)
  names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null",
               "Likelihood.ratio.test", "Number.of.observations",
               "Messages", "Warnings")
  return(V)

}

extract.nagelkerke.r2 <- function(x) {
  # input nagelkerke() output, extracts and outputs relevant r2

  r2 <- x$Pseudo.R.squared.for.model.vs.null[3, 1]
  r2_percentage_points <- r2 * 100

  return(r2_percentage_points)
}


########################################################
###################### Data paths ######################
########################################################

# scores
score_wd <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/public/scores/"
score_files <- list.files(score_wd)
HRS_files_all <- score_files[grep("HRS", score_files)]
score_names <- gsub("PGS_HRS2_", "", gsub("_LDpred_p1.txt", "", HRS_files_all))
score_names <- gsub("-.*", "", score_names)


# phenos
pheno_wd <- "input/HRS"
pheno_files <- list.files(pheno_wd)
pheno_names <- gsub(".pheno", "", pheno_files)

# overlapping names
both_names <- pheno_names[pheno_names %in% score_names]

# check for correct pheno name overlap:
score_names[!(score_names %in% both_names)]; both_names; pheno_names

# continue overlapping names, to get loop working
pheno_names <- both_names

# scores-phenos crosswalk
score_pheno_crosswalk_path <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.dta"
score_pheno_crosswalk_data <- read.dta(score_pheno_crosswalk_path) %>%
  mutate(IID = as.numeric(LOCAL_ID),
         HHID = as.numeric(HHID),
         PN = as.numeric(PN))

# PCs and related crosswalk
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
  inc_r2_public = numeric()
)

# options for prediction loop
iterations <- 1000
binary_phenotypes <- c(
  "ADHD",
  "ASTECZRHI",
  "ASTHMA",
  "EVERSMOKE",
  "HAYFEVER",
  "MIGRAINE",
  "NEARSIGHTED"
)
use_nagelkerke <- FALSE

for (i in 1:length(pheno_names)){

  local_df <- data.frame(phenotype = character(1))

  pheno <- pheno_names[i]
  local_df$phenotype <- pheno

  pheno_path <- paste0(getwd(), "/input/HRS/", pheno, ".pheno")
  local_score_path <- HRS_files_all[grep(pheno, HRS_files_all,
                                           fixed=TRUE)]

  score_path <- paste0(score_wd, local_score_path)

  pheno_data <- fread(pheno_path)
  score_data <- fread(score_path) %>%
    rename(score = contains(substr(pheno, 1, 3)))

  if ("hhidpn" %in% colnames(pheno_data)){ # if needing to join on hhidpn variables
    score_pheno_crosswalk_data <- score_pheno_crosswalk_data %>% mutate(hhidpn = (1000 * HHID) + PN)
    data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
      inner_join(pheno_data, by="hhidpn") %>%
      inner_join(PCs_data, by="IID") %>%
      drop_na(phenotype, score)
  } else { # else, if joining on hhid and pn
      if ("hhid" %in% colnames(pheno_data)) {
        pheno_data <- pheno_data %>% mutate(HHID = hhid, PN = pn)
      }
      data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
        inner_join(pheno_data, by=c("HHID", "PN")) %>%
        inner_join(PCs_data, by="IID") %>%
        drop_na(phenotype, score)
    }

  # get PCs
  PCs <- str_c(paste0("pc", 1:10), collapse=" + ")

  # regress with and without scores, extract r^2
  formula_no_score <- as.formula(paste0("phenotype ~ ", PCs))
  formula_with_score <- as.formula(paste0("phenotype ~ score + ", PCs))
  reg_no_score <- lm(formula_no_score, data)
  reg_with_score <- lm(formula_with_score, data)
  if (pheno %in% binary_phenotypes & use_nagelkerke) {
    # if outcome binary, calculate nagelkerke r^2

    r2_no_score <- nagelkerke(resampled_no_score_reg_fit)
    Sys.sleep(0.1)
    r2_with_score <- nagelkerke(resampled_with_score_reg_fit)

    local_df$r2_no_score <- extract.nagelkerke.r2(r2_no_score)
    Sys.sleep(0.1)
    local_df$r2_with_score <- extract.nagelkerke.r2(r2_with_score)
  } else {
    # if outcome not binary, calculate standard r^2

    local_df$r2_no_score <- summary(reg_no_score)$r.squared * 100
    local_df$r2_with_score <- summary(reg_with_score)$r.squared * 100
  }
  local_df$r2_inc <- local_df$r2_with_score - local_df$r2_no_score

  # bootstrap
  r2_inc_list <- c()

  # iterate for bootstrapping
  for (j in 1:iterations) {

    # print iteration and phenotype
    if (j %% 100 == 0){
      print(paste0("Iteration: ", j, ". Phenotype: ", pheno, "."))
    }

    # resample, regress, add r2
    resampled_data <- sample_n(data, nrow(data), replace=TRUE)
    resampled_no_score_reg_fit <- lm(formula_no_score, resampled_data)
    resampled_with_score_reg_fit <- lm(formula_with_score, resampled_data)

    if (pheno %in% binary_phenotypes & use_nagelkerke) {
      # if outcome binary, calculate nagelkerke r^2

      if (j == 1){
        print(paste0("Using Nagelkerke R^2 for phenotype: ", pheno, "."))
      }

      r2_no_score <- nagelkerke(resampled_no_score_reg_fit)
      Sys.sleep(0.1)
      r2_with_score <- nagelkerke(resampled_with_score_reg_fit)

      r2_no_score <- extract.nagelkerke.r2(r2_no_score)
      Sys.sleep(0.1)
      r2_with_score <- extract.nagelkerke.r2(r2_with_score)
    } else {
      # if outcome not binary, calculate standard r^2
      resampled_no_score_reg <- summary(resampled_no_score_reg_fit)
      resampled_with_score_reg <- summary(resampled_with_score_reg_fit)

      r2_no_score <- resampled_no_score_reg$r.squared * 100
      r2_with_score <- resampled_with_score_reg$r.squared * 100
    }
    r2_inc_list[j] <- r2_with_score - r2_no_score

  }

  # save bootstrapped r2s
  local_bootstrap_df <- data.frame(
    iteration = 1:1000,
    phenotype = pheno
  )
  local_bootstrap_df$inc_r2_public <- r2_inc_list
  bootstrap_df <- rbind(bootstrap_df, local_bootstrap_df)

  # generate quantities of interest
  r2_inc_low  <- quantile(r2_inc_list, .025)[[1]]
  r2_inc_high <- quantile(r2_inc_list, .975)[[1]]
  r2_inc_mean <- mean(r2_inc_list)

  # print r2 confidence interval
  print(paste0("Incremental R^2 (95% CI): ", r2_inc_mean, " [", r2_inc_low, " : ", r2_inc_high, "]"))

  ### bootstrap function
  local_df$r2_inc_lower <- r2_inc_low
  local_df$r2_inc_upper <- r2_inc_high

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
    "/output/HRS_bootstraps_public.txt"
  )
)

# overwrite main file
fwrite(
  df,
  paste0(
    getwd(),
    "/output/HRS_phenotypes_r2_public.txt"
  )
)
