#!/bin/bash

source paths11

cd $mainDir/derived_data/11_Figures/

Rscript $mainDir/code/11_Figures/11.2.0_gendata_supptbl.R $mainDir
Rscript $mainDir/code/11_Figures/11.2.1_plotting5datasets_single.R $mainDir
Rscript $mainDir/code/11_Figures/11.2.2_plotting5datasets_multi.R $mainDir
Rscript $mainDir/code/11_Figures/11.2.3_plotting5datasets_diff.R $mainDir