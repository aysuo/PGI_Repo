source $PGI_Repo/code/paths
source $PGI_Repo/code/10_Prediction/paths_PGIrepo_withinFam

Rscript 11.1_plotting_IncR2.R \
    $PGI_Repo/derived_data/10_Prediction/UKB3_EUR_r2.txt \
    $PGI_Repo/derived_data/10_Prediction/HRS_EUR_r2.txt \
    $PGI_Repo/derived_data/10_Prediction/WLS_EUR_r2.txt \
    $PGI_RepoV2/doc/PhenoList/fullPheno_name_list.txt \
    $PGI_Repo/derived_data/11_Figures/output \
    NA
#------

sh 11.2_Fig3_causal_vs_population_effects_UKB.R

sh 11.3_FigS3_ratio_causal_vs_pop_effects_UKB.R

sh 11.4_FigS4_causal_vs_population_effects_WLS.R

sh 11.5_FigS5_ratio_causal_vs_pop_effects_WLS.R

sh 11.6_Fig4_causalVsPopEffects_byDemography_UKB.R

sh 11.7_FigS6_IncrR2_ratio_byDemography_UKB.R