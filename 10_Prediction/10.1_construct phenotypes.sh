#!/bin/bash

# WLS
stata -b 10.1.0_save_WLS.do
Rscript 10.1.1_construct_WLS_phenotypes.R 


# HRS
Rscript 10.1.2_construct_HRSRAND_phenotypes.R
stata -b 10.1.3_construct_HRS_SWB_part1.do
Rscript 10.1.4_construct_HRS_SWB_part2.R
stata -b 10.1.5_construct_HRS_CIDI.do
