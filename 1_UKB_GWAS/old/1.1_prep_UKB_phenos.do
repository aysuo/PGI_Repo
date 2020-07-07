clear all
cd "/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS"
local crosswalk="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_imp_v2.ukb1142_imp_chr1_v2_s487398.crosswalk"
local partition_data="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt"
local pheno_data_1="/disk/genetics2/ukb/orig/UKBv2/pheno/9690_update/ukb9690_update_jun14_2018.dta"
local covar_data="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_sqc_v2_combined_header.dta"
local pheno_data_2="/disk/genetics2/ukb/orig/UKBv2/pheno/22392/ukb22392.dta"

log using "/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/1_UKB_GWAS/1.1_prep_UKB_phenos.do.log", replace
display "$S_DATE $S_TIME"



******************************************************
**                 PREPARE PHENOTYPES UKB           **
**                  PGS Repo project                **
**                                                  **
**                Casper Burik 12 Feb 2018          **
**              Last edit Casper: 22 Feb 2020       **
**              Last edit Aysu: 3 March 2020        **
******************************************************

** To do: Risk, Alzheimer's proxy
** Check which phenotypes updated: migraine, hayfever, depression
** Get phenotypes: AUDIT, Cannabis, Age first menses, Eczema?, n_22040 // MET minutes per week
** Get phenotypes: 22126 (for hayfever) 22127 (for asthma)

set more off
set maxvar 32000


use `pheno_data_1'

keep n_50_* ///  Height 
    n_23104_* /// BMI
    n_20127_0_0 /// Neuroticism scoreplot
    n_1180_* /// Morning person
    n_6160* /// Religious group activity
    n_4526*  /// Happiness / Subjective Wellbeing
    n_20160* /// Ever smoked
    n_2178* /// Self-rated Health
    n_2754* /// Age at First birth
    n_2734* /// Number of children (women)
    n_2405* /// Number of children (men)
    n_6152* /// Other diagnoses (Hayfever)
    s_41202_* s_41204_* /// ICD10
    n_20002_* /// Non-cancer illness
    /// n_3872* /// Age at First birth single child
    n_4559_* /// Family satisfaction
    n_4570_* /// Friendship satisfaction
    n_4537_* /// Work satisfaction
    n_4581_* /// Financial satisfaction
    n_2020_* /// loneliness
    n_3456_* n_2887_* n_6183_* /// Cigarettes per day
    n_19?0_* ///
    n_2000_* ///
    n_2010_* ///
    n_2020_* ///
    n_2030_* ///
    n_1558_* /// drinking frequency
    n_1568_* /// red wine per week
    n_1578_* /// white wine and champagne
    n_1588_* /// beer and cider per week
    n_1598_* /// spirits per week
    n_1608_* /// fortified wine per week
    n_5364_* /// other alcoholic drinks per week LESS PEOPLE
    n_4407_* /// red wine per month
    n_4418_* /// white and champagne
    n_4429_* /// beer
    n_4440_* /// spirits
    n_4451_* /// fortified wine
    n_4462_* /// other
    n_20117_* /// EVERDRINK status (derived)
    n_1239_* /// current smoking
    n_1249_* /// past smoking
    n_2644_* /// light smoking
    n_21000_* /// // Ethnicity
    n_34_0_0 n_22200_0_0 /// birth year
    n_22010_* /// // Bad genotype n_22027_0_0
    n_eid // ID

gen IID=n_eid
merge 1:1 IID using `covar_data', nogen keep(match) 
merge 1:1 n_eid using `pheno_data_2', nogen keep(match)

*** DROP OBSERVATIONS ***
// Withdrawn individuals:
foreach id in 1038534 1048252 1085993 1136386 1138937 1142637 1154859 1262850 1283111 1283616 ///
    1300553 1364378 1370374 1371898 1378890 1405916 1425122 1431973 1437530 1451103 1481466 ///
    1533928 1535412 1548599 1560258 1560429 1604306 1626372 1750924 1799527 1872177 1890503 ///
    1977732 2002238 2040924 2099267 2119447 2165727 2238632 2242976 2323750 2332131 2472313 ///
    2474972 2512665 2518039 2586585 2679822 2687632 2702195 2723486 2784782 2791953 2793515 ///
    2796360 2947991 2994818 3091905 3098976 3157343 3162399 3178848 3199605 3258519 3314334 ///
    3373335 3393798 3403770 3447432 3466177 3553318 3563210 3580955 3590961 3623224 3639236 ///
    3641934 3656893 3672665 3675369 3739874 3789586 3877018 3880771 3931774 3939246 3962080 ///
    3964248 3975476 4012699 4025282 4070440 4084812 4155673 4170030 4183522 4337448 4405797 ///
    4419037 4598042 4693228 4701857 4732737 4819063 4823640 4884464 4893886 4962254 4980287 ///
    4992864 5070668 5087142 5106856 5137728 5205453 5211729 5225718 5266776 5267820 5291977 ///
    5305498 5312601 5344717 5356287 5368810 5369682 5388606 5401275 5416501 5417891 5464094 ///
    5487522 5554364 5569334 5571224 5654288 5689279 5697575 5702369 5730419 5750297 5762972 ///
    5770853 5801079 5878770 5914593 5950775 5958708 6008079{
        drop if n_eid == `id'
    }


*** QC ***
ren Sex SEX

* Undefined sex
drop if SEX == 0 
drop if missing(SEX)

ren Batch BATCH
drop if missing(BATCH)
drop if InferredGender != SubmittedGender

* Bad genotypes, initial release
drop if n_22010_0_0 == 1 
* Bad genotypes, second release
drop if hetmissingoutliers == 1
drop if missing(PC1-PC40)
drop if n_eid < 0 
drop if putativesexchromosomeaneuploidy == 1

* Non-Europeans
drop if PC1 > 0 

*** Create ETHN  ***
gen ETHN =  n_21000_0_0
    replace ETHN = n_21000_1_0 if missing(ETHN)
    replace ETHN = n_21000_2_0 if missing(ETHN)
    drop if missing(ETHN)
    ** Keep White (value = 1) / British (value = 1001) / Irish (value = 1002) / Any other white background (value = 1003)
    drop if (ETHN != 1 & ETHN != 1001 & ETHN != 1002 & ETHN != 1003)

*** Create BYEAR ***
gen BYEAR = (n_34_0_0 - 1900)/10
    replace BYEAR = (n_22200_0_0 - 1900)/10 if missing(BYEAR)
    drop if missing(BYEAR)



**********************************************************************************
******************************* GENERATE VARS ************************************
**********************************************************************************

**********************************************************
************************ HEIGHT **************************
**********************************************************
egen HEIGHT =  rmean(n_50_*)
**********************************************************


**********************************************************
************************* BMI ****************************
**********************************************************
egen BMI = rmean(n_23104_*)
**********************************************************


**********************************************************
********************* NEUROTICISM ************************
**********************************************************
ren n_20127_0_0 NEURO
**********************************************************


**********************************************************
********************* MORNING PERSON *********************
**********************************************************
* scale 1-4, was reverse coded
gen MORNING_0 =  5 - n_1180_0_0 if n_1180_0_0 > 0
* <0 means missing 
gen MORNING_1 = 5 - n_1180_1_0 if n_1180_1_0 > 0
gen MORNING_2 = 5 - n_1180_2_0 if n_1180_2_0 > 0
egen MORNING = rmean(MORNING_*)
**********************************************************


**********************************************************
*********************** RELIGIOSITY **********************
**********************************************************
* Religious = 1 if answered 'Religious group' at least once

gen RELIG = 0 if (n_6160_0_0 != -3 & !missing(n_6160_0_0))
replace RELIG = 0 if (missing(RELIG) & !missing(n_6160_1_0) & n_6160_1_0 != -3)
replace RELIG = 0 if (missing(RELIG) & !missing(n_6160_2_0) & n_6160_2_0 != -3)

forval i = 0/2 {
    forval j = 0/4 {
        replace RELIG = 1 if n_6160_`i'_`j' == 3
    }
}
**********************************************************


**********************************************************
*************************** SWB **************************
**********************************************************
* scale 1-6, was reverse coded
gen SWB_0 = 7 - n_4526_0_0 if n_4526_0_0 > 0 
gen SWB_1 = 7 - n_4526_1_0 if n_4526_1_0 > 0
gen SWB_2 = 7 - n_4526_2_0 if n_4526_2_0 > 0
egen SWB = rmean(SWB_*)
**********************************************************

*rerun
**********************************************************
*********************** DEPRESSION ***********************
**********************************************************
* replace 'do not knows'
forval i=2/9 {
    replace n_19`i'0_0_0 = . if n_19`i'0_0_0 < 0
    replace n_19`i'0_1_0 = . if n_19`i'0_1_0 < 0
    replace n_19`i'0_2_0 = . if n_19`i'0_2_0 < 0
}

forval i=0/3 {
    replace n_20`i'0_0_0 = . if n_20`i'0_0_0 < 0
    replace n_20`i'0_1_0 = . if n_20`i'0_1_0 < 0
    replace n_20`i'0_2_0 = . if n_20`i'0_2_0 < 0
}

* take means of repeat observations
forval i=2/9 {
    egen n_19`i'0 = rmean(n_19`i'0_*)
}
forval i=0/3 {
    egen n_20`i'0 = rmean(n_20`i'0_*)
}

* sum
gen DEPRESS = n_1920 + n_1930 + n_1940 + n_1950 + n_1960 + n_1970 + n_1980 ///
            + n_1990 + n_2000 + n_2010 + n_2020 + n_2030
**********************************************************


**********************************************************
****************** FAMILY SATISFACTION *******************
**********************************************************
gen FAM_SATISF_0 = 7 - n_4559_0_0 if n_4559_0_0 > 0
gen FAM_SATISF_1 = 7 - n_4559_1_0 if n_4559_1_0 > 0
gen FAM_SATISF_2 = 7 - n_4559_2_0 if n_4559_2_0 > 0
egen FAM_SATISF = rmean(FAM_SATISF_*)
**********************************************************


**********************************************************
****************** FRIEND SATISFACTION *******************
**********************************************************
gen FRIEND_SATISF_0 = 7 - n_4570_0_0 if n_4570_0_0 > 0
gen FRIEND_SATISF_1 = 7 - n_4570_1_0 if n_4570_1_0 > 0
gen FRIEND_SATISF_2 = 7 - n_4570_2_0 if n_4570_2_0 > 0
egen FRIEND_SATISF = rmean(FRIEND_SATISF_*)
**********************************************************


**********************************************************
***************** FINANCIAL SATISFACTION *****************
**********************************************************
gen FIN_SATISF_0 = 7 - n_4581_0_0 if n_4581_0_0 > 0
gen FIN_SATISF_1 = 7 - n_4581_1_0 if n_4581_1_0 > 0
gen FIN_SATISF_2 = 7 - n_4581_2_0 if n_4581_2_0 > 0
egen FIN_SATISF = rmean(FIN_SATISF_*)
**********************************************************


**********************************************************
******************* WORK SATISFACTION ********************
**********************************************************
gen WORK_SATISF_0 = 7 - n_4537_0_0 if (n_4537_0_0 > 0 & n_4537_0_0 != 7)
* 7 is unemployed
gen WORK_SATISF_1 = 7 - n_4537_1_0 if (n_4537_1_0 > 0 & n_4537_1_0 != 7)
gen WORK_SATISF_2 = 7 - n_4537_2_0 if (n_4537_2_0 > 0 & n_4537_2_0 != 7)
egen WORK_SATISF = rmean(WORK_SATISF*)
**********************************************************


**********************************************************
********************** LONELINESS ************************
**********************************************************
forval i = 0/2 {
    replace n_2020_`i'_0=. if n_2020_`i'_0<0
}

egen LONELY = rmean(n_2020_0_0 n_2020_1_0 n_2020_2_0) 
**********************************************************


**********************************************************
*********************** SELFHEALTH ***********************
**********************************************************
* scale 1-4, was reverse coded
gen HEALTH_0 = 5 - n_2178_0_0 if n_2178_0_0 > 0 
gen HEALTH_1 = 5 - n_2178_1_0 if n_2178_1_0 > 0
gen HEALTH_2 = 5 - n_2178_2_0 if n_2178_2_0 > 0
egen HEALTH = rmean(HEALTH_*)
**********************************************************


**********************************************************
************************** ASTHMA ************************
**********************************************************
** Update when we have data field 22127!!

gen ASTHMA = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
    (!missing(n_6152_2_0)) &  n_6152_2_0!=-3) 


forval i = 0/379 {
    replace ASTHMA = 1 if (s_41202_0_`i' == "J450" | s_41202_0_`i' == "J451" | ///
            s_41202_0_`i' == "J458" |  s_41202_0_`i' == "J459" | ///
            s_41202_0_`i' == "J46")
}

forval i = 0/434 {
    replace ASTHMA = 1 if (s_41204_0_`i' == "J450" | s_41204_0_`i' == "J451" | ///
            s_41204_0_`i' == "J458" |  s_41204_0_`i' == "J459"   | ///
            s_41204_0_`i' == "J46")
}

forval i = 0/28 {
    replace ASTHMA = 1 if (n_20002_0_`i' == 1111 )
}

forval i = 0/15 {
    replace ASTHMA = 1 if (n_20002_1_`i' == 1111 )
}

forval i = 0/16 {
    replace ASTHMA = 1 if (n_20002_2_`i' == 1111 )
}

forval i = 0/3 {
    forval j = 0/2 {
        replace ASTHMA = 1 if n_6152_`j'_`i' == 8
    }
}
replace ASTHMA = 1 if n_6152_0_4 == 8



**********************************************************


*rerun
**********************************************************
*********************** HAYFEVER *************************
**********************************************************
** Update when we get data field 22126

gen HAYFVR = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3 & n_6152_0_0!=9) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3 & n_6152_0_0!=9) | ///
    (!missing(n_6152_2_0) & n_6152_2_0!=-3 & n_6152_0_0!=9)) 


forval i = 0/379 {
    replace HAYFVR = 1 if (s_41202_0_`i' == "J301" | s_41202_0_`i' == "J302" | ///
        s_41202_0_`i' == "J303" |  s_41202_0_`i' == "J304")
}

forval i = 0/434 {
    replace HAYFVR = 1 if (s_41204_0_`i' == "J301" | s_41204_0_`i' == "J302" | ///
        s_41204_0_`i' == "J303" |  s_41204_0_`i' == "J304")
}

forval i = 0/28 {
    replace HAYFVR = 1 if (n_20002_0_`i' == 1387 )
}

forval i = 0/15 {
    replace HAYFVR = 1 if (n_20002_1_`i' == 1387 )
}

forval i = 0/16 {
    replace HAYFVR = 1 if (n_20002_2_`i' == 1387 )
}
**********************************************************



**********************************************************
************************ ECZEMA **************************
**********************************************************
gen ECZEMA = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3 & n_6152_0_0!=9) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3 & n_6152_0_0!=9) | ///
    (!missing(n_6152_2_0) & n_6152_2_0!=-3 & n_6152_0_0!=9)) 


forval i = 0/379 {
    replace ECZEMA = 1 if (s_41202_0_`i' == "L20" | s_41202_0_`i' == "L208" | ///
        s_41202_0_`i' == "J209")
}

forval i = 0/434 {
    replace ECZEMA = 1 if (s_41204_0_`i' == "L20" | s_41204_0_`i' == "L208" | ///
        s_41204_0_`i' == "J209")
}

forval i = 0/28 {
    replace ECZEMA = 1 if (n_20002_0_`i' == 1452 )
}

forval i = 0/15 {
    replace ECZEMA = 1 if (n_20002_1_`i' == 1452 )
}

forval i = 0/16 {
    replace ECZEMA = 1 if (n_20002_2_`i' == 1452 )
}
**********************************************************


**********************************************************
******************* HAYFEVER & ECZEMA ********************
**********************************************************
gen ECZRHI = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
    (!missing(n_6152_2_0) &  n_6152_2_0!=-3)) 


forval i = 0/2 {
    forval j = 0/3 {
        replace ECZRHI = 1 if n_6152_`i'_`j' == 9
    }
}
    
replace ECZRHI = 1 if n_6152_0_4 == 9

forval i = 0/379 {
    replace ECZRHI = 1 if (s_41202_0_`i' == "J301" | s_41202_0_`i' == "J302" | ///
        s_41202_0_`i' == "J303" |  s_41202_0_`i' == "J304" | s_41202_0_`i' == "L20" | /// 
        s_41202_0_`i' == "L208" |  s_41202_0_`i' == "J209")
}

forval i = 0/434 {
    replace ECZRHI = 1 if (s_41204_0_`i' == "J301" | s_41204_0_`i' == "J302" | ///
        s_41204_0_`i' == "J303" |  s_41204_0_`i' == "J304" | s_41204_0_`i' == "L20" | ///
        s_41204_0_`i' == "L208" |  s_41204_0_`i' == "J209")
}

forval i = 0/28 {
    replace ECZRHI = 1 if (n_20002_0_`i' == 1387 | n_20002_0_`i' == 1452)
}

forval i = 0/15 {
    replace ECZRHI = 1 if (n_20002_1_`i' == 1387 | n_20002_1_`i' == 1452)
}

forval i = 0/16 {
    replace ECZRHI = 1 if (n_20002_2_`i' == 1387 | n_20002_2_`i' == 1452)
}
**********************************************************



**********************************************************
************** ASTHMA / HAYFEVER/ECZEMA ******************
**********************************************************
gen ASTECZRHI=.
replace ASTECZRHI=1 if (ASTHMA==1 | ECZRHI==1)
replace ASTECZRHI=0 if (ASTHMA==0 & ECZRHI==0)

**********************************************************



**********************************************************
*********************** EVERSMOKE ************************
**********************************************************

forval i = 0/2 {
    gen current_`i' = 1 if n_1239_`i'_0 == 1 
    replace current_`i' = 0 if n_1239_`i'_0 == 0 |  n_1239_`i'_0 == 2

    gen past_`i' = 1 if n_1249_`i'_0 == 1 | (n_1249_`i'_0 == 2 & n_2644_`i'_0 == 1) | (n_1249_`i'_0 == 3 & n_2644_`i'_0 == 1)
    replace past_`i' = 0 if n_1249_`i'_0 == 4 | (n_1249_`i'_0 == 2 & n_2644_`i'_0 == 0) | (n_1249_`i'_0 == 3 & n_2644_`i'_0 == 0)

    egen EVERSMOKE_`i'=rmax(current_`i' past_`i')
}

egen EVERSMOKE = rmax(EVERSMOKE_*)

**********************************************************


*rerun
**********************************************************
******************* CIGARETTES PER DAY *******************
**********************************************************

forval i = 0/2 {
    replace n_3456_`i'_0=. if n_3456_`i'_0 < 1 | n_3456_`i'_0 > 100
    replace n_2887_`i'_0=. if n_2887_`i'_0 < 1 | n_2887_`i'_0 > 100
    replace n_6183_`i'_0=. if n_6183_`i'_0 < 1 | n_6183_`i'_0 > 100
}

egen CPD = rmean(n_3456_* n_2887_* n_6183_*)

* GSCAN binning:
* a. 1 = 1-5
* b. 2 = 6-15
* c. 3 = 16-25
* d. 4 = 26-35
* e. 5 = 36+

gen CPDbins=1 if CPD > 0 & CPD <= 5
replace CPDbins=2 if CPD > 5 & CPD <= 15
replace CPDbins=3 if CPD > 15 & CPD <= 25
replace CPDbins=4 if CPD > 25 & CPD <= 35
replace CPDbins=5 if CPD > 35 & CPD != .
**********************************************************


*rerun
**********************************************************
******************** DRINKS PER WEEK *********************
**********************************************************
forval i = 0/2 {
    forval j=56/60 {
        replace n_1`j'8_`i'_0 = . if n_1`j'8_`i'_0 < 0
    }

    replace n_5364_`i'_0 = . if n_5364_`i'_0 < 0
    replace n_4407_`i'_0 = . if n_4407_`i'_0 < 0
    replace n_4418_`i'_0 = . if n_4418_`i'_0 < 0
    replace n_4429_`i'_0 = . if n_4429_`i'_0 < 0
    replace n_4440_`i'_0 = . if n_4440_`i'_0 < 0
    replace n_4451_`i'_0 = . if n_4451_`i'_0 < 0
    replace n_4462_`i'_0 = . if n_4462_`i'_0 < 0
}

forval i = 0/2 {
    * Set to sum of drinks if not less frequent than once or twice per week 
    egen AW_`i' = rowtotal(n_1568_`i'_0 n_1578_`i'_0 n_1588_`i'_0 n_1598_`i'_0 n_1608_`i'_0 n_5364_`i'_0) if  (n_1568_`i'_0!=. | n_1578_`i'_0!=. | n_1588_`i'_0!=. | n_1598_`i'_0!=. |  n_1608_`i'_0!=. |  n_5364_`i'_0!=.)

 
    * If less frequent than once or twice per week, replace with rescaled monthly value 
    egen AM_`i' = rowtotal(n_4407_`i'_0 n_4418_`i'_0 n_4429_`i'_0 n_4440_`i'_0 n_4451_`i'_0 n_4462_`i'_0) if (n_4407_`i'_0!=. | n_4418_`i'_0!=. | n_4429_`i'_0!=. | n_4440_`i'_0!=. | n_4451_`i'_0!=. | n_4462_`i'_0!=.)

    replace AW_`i'= AM_`i' / 4 if (n_1558_`i'_0 == 4 | n_1558_`i'_0 == 5 )
    
    * Set to 0 if never drinks 
    replace AW_`i'= 0 if n_1558_`i'_0 == 6
}

egen DPW = rmean(AW_*)

**********************************************************



**********************************************************
************************** COPD **************************
**********************************************************
* ICD10 codes received from here:
* https://www.cigna.com/static/docs/starplus/icd10-copd.pdf

gen COPD = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
    (!missing(n_6152_2_0)) &  n_6152_2_0!=-3) 

forval i = 0/379 {
    replace COPD = 1 if (s_41202_0_`i' == "J440" | s_41202_0_`i' == "J441" | ///
             s_41202_0_`i' == "J448" |  s_41202_0_`i' == "J449" | s_41202_0_`i' == "J410" | ///
             s_41202_0_`i' == "J411" | s_41202_0_`i' == "J418" | s_41202_0_`i' == "J42" | ///
             s_41202_0_`i' == "J430" | s_41202_0_`i' == "J431" | s_41202_0_`i' == "J432" | ///
             s_41202_0_`i' == "J438" | s_41202_0_`i' == "J439" | s_41202_0_`i' == "J982" | ///
             s_41202_0_`i' == "J983")
}

forval i = 0/434 {
    replace COPD = 1 if (s_41204_0_`i' == "J440" | s_41204_0_`i' == "J441" | ///
             s_41204_0_`i' == "J448" |  s_41204_0_`i' == "J449" | s_41204_0_`i' == "J410" | ///
             s_41204_0_`i' == "J411" | s_41204_0_`i' == "J418" | s_41204_0_`i' == "J42" | ///
             s_41204_0_`i' == "J430" | s_41204_0_`i' == "J431" | s_41204_0_`i' == "J432" | ///
             s_41204_0_`i' == "J438" | s_41204_0_`i' == "J439" | s_41204_0_`i' == "J982" | ///
             s_41204_0_`i' == "J983")
}

forval i = 0/28 {
    replace COPD = 1 if (n_20002_0_`i' == 1112 | n_20002_0_`i' == 1113)
}

forval i = 0/15 {
    replace COPD = 1 if (n_20002_1_`i' == 1112 | n_20002_1_`i' == 1113)
}

forval i = 0/16 {
    replace COPD = 1 if (n_20002_2_`i' == 1112 | n_20002_2_`i' == 1113)
}

forval i = 0/3 {
    forval j = 0/2 {
        replace COPD = 1 if n_6152_`j'_`i' == 6
    }
}
replace COPD = 1 if n_6152_0_4 == 6

**********************************************************


*rerun
**********************************************************
************************* MIGRAINE ***********************
**********************************************************
gen MIGRAINE = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
        !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0))


forval i = 0/379 {
    replace MIGRAINE = 1 if (s_41202_0_`i' == "G430" | s_41202_0_`i' == "G431" | ///
            s_41202_0_`i' == "G432" |  s_41202_0_`i' == "G433" | ///
            s_41202_0_`i' == "G438" | s_41202_0_`i' == "G439")
}

forval i = 0/434 {
    replace MIGRAINE = 1 if (s_41204_0_`i' == "G430" | s_41204_0_`i' == "G431" | ///
            s_41204_0_`i' == "G432" |  s_41204_0_`i' == "G433"   | ///
            s_41204_0_`i' == "G438" | s_41204_0_`i' == "G439")
}

forval i = 0/28 {
    replace MIGRAINE = 1 if (n_20002_0_`i' == 1265 )
}

forval i = 0/15 {
    replace MIGRAINE = 1 if (n_20002_1_`i' == 1265 )
}

forval i = 0/16 {
    replace MIGRAINE  = 1 if (n_20002_2_`i' == 1265 )
}

**********************************************************



**********************************************************
******************** AGE AT FIRST BIRTH ******************
**********************************************************
* minimum age of all times answered
gen AFB_0 = n_2754_0_0 if n_2754_0_0 > 0
gen AFB_1 = n_2754_1_0 if n_2754_1_0 > 0
gen AFB_2 = n_2754_2_0 if n_2754_2_0 > 0
gen AFB_0b = n_3872_0_0 if n_3872_0_0 > 0
gen AFB_1b = n_3872_1_0 if n_3872_1_0 > 0
gen AFB_2b = n_3872_2_0 if n_3872_2_0 > 0
egen AFB = rmin(AFB_*)
**********************************************************



**********************************************************
******************** NUMBER EVER BORN ********************
**********************************************************
* max of number all reported children
gen N_EVERBORN_WOMEN_0 = n_2734_0_0 if n_2734_0_0 >= 0
gen N_EVERBORN_WOMEN_1 = n_2734_1_0 if n_2734_1_0 >= 0
gen N_EVERBORN_WOMEN_2 = n_2734_2_0 if n_2734_2_0 >= 0
egen N_EVERBORN_WOMEN = rmax(N_EVERBORN_WOMEN_*)

gen N_EVERBORN_MEN_0 = n_2405_0_0 if n_2405_0_0 >= 0
gen N_EVERBORN_MEN_1 = n_2405_1_0 if n_2405_1_0 >= 0
gen N_EVERBORN_MEN_2 = n_2405_2_0 if n_2405_2_0 >= 0
egen N_EVERBORN_MEN = rmax(N_EVERBORN_MEN_*)


**********************************************************************
**********************************************************************

*** SAVE FULL DATASET ***
keep n_eid SEX BATCH BYEAR HEIGHT BMI NEURO MORNING RELIG SWB HEALTH AFB N_EVERBORN_MEN  ///
    N_EVERBORN_WOMEN HAYFVR ECZEMA ECZRHI ASTECZRHI FAM_SATISF FRIEND_SATISF WORK_SATISF FIN_SATISF ///
    LONELY CPD CPDbins COPD ASTHMA MIGRAINE DEPRESS DPW EVERSMOKE PC*

save "tmp/pgs_repo.dta", replace


**********************************************************************************
**********************************************************************************

****************************************
**** Split into 3 equal partitions *****
****************************************
clear
import delimited `partition_data', delim("\t")
merge 1:1 n_eid using "tmp/pgs_repo.dta", keep(match using)

* Get maximum order in IDs_assignPartition_ordered.txt after merging 
qui sum partition_order

* For observations not in brain sample + relateds list, assign random value greater than maximum order to partition_order 
local maxorder =`r(max)'+1
replace partition_order=runiform()+`maxorder' if partition_order==.

* Assign into partitions based on partition order
sort partition_order
gen partition = ceil(3 * _n/_N)
save "tmp/pgs_repo.dta", replace

****************************************
********* Convert to GWAS IDs **********
****************************************

clear 
import delimited `crosswalk', delim(" ")
ren (v1 v2) (IID n_eid)
destring n_eid, replace ignore(" ")
merge 1:1 n_eid using "tmp/pgs_repo.dta", nogen
gen FID=IID


****************************************
********* RESIDUALIZE & EXPORT *********
****************************************

gen BYEAR2 = BYEAR*BYEAR
gen BYEAR3 = BYEAR2*BYEAR
gen SEXxBYEAR = SEX*BYEAR
gen SEXxBYEAR2 = SEX*BYEAR2
gen SEXxBYEAR3 = SEX*BYEAR3


foreach partition in 1 2 3{
    foreach var of varlist FAM_SATISF FRIEND_SATISF WORK_SATISF FIN_SATISF ///
        LONELY CPD CPDbins COPD ASTHMA MIGRAINE DEPRESS DPW EVERSMOKE ///
        HEIGHT BMI NEURO MORNING RELIG SWB HEALTH HAYFVR ECZEMA ECZRHI ASTECZRHI{
        
        qui xi:reg `var' BYEAR* SEX* i.BATCH PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }


    * sex specific phenotypes:
    foreach var of varlist AFB N_EVERBORN_MEN N_EVERBORN_WOMEN {
        
        qui xi:reg `var' BYEAR* i.BATCH PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }
}


*** END ***
log close
