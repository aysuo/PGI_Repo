#!/bin/bash

cd $PGI_Repo/derived_data/11_Figures/

Rscript $PGI_Repo/code/11_Figures/11.2.0_gendata_supptbl.R $PGI_Repo
Rscript $PGI_Repo/code/11_Figures/11.2.1_plotting5datasets_single.R $PGI_Repo
Rscript $PGI_Repo/code/11_Figures/11.2.2_plotting5datasets_multi.R $PGI_Repo
Rscript $PGI_Repo/code/11_Figures/11.2.3_plotting5datasets_diff.R $PGI_Repo