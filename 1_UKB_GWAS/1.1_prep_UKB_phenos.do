clear all
cd "/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS"
local crosswalk="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_imp_v2.ukb1142_imp_chr1_v2_s487398.crosswalk"
local partition_data="/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt"
local pheno_data_1="/disk/genetics2/ukb/orig/UKBv2/pheno/9690_update/ukb9690_update_jun14_2018.dta"
local covar_data="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_sqc_v2_combined_header.dta"
local pheno_data_2="/disk/genetics2/ukb/orig/UKBv2/pheno/22392/ukb22392.dta"
local pheno_data_3="/disk/genetics/ukb/data/41490_alc_cannabis_menarche_myopia/ukb41490.dta"

log using "/disk/genetics4/PGS/Aysu/PGS_Repo_pipeline/code/1_UKB_GWAS/1.1_prep_UKB_phenos.do.log", replace
display "$S_DATE $S_TIME"



******************************************************
**                 PREPARE PHENOTYPES UKB           **
**                  PGS Repo project                **
**                                                  **
**                Casper Burik 12 Feb 2018          **
**              Last edit Casper: 22 Feb 2020       **
**              Last edit Aysu: 14 April 2020        **
******************************************************

** Get phenotypes: Eczema?, 22126 (for hayfever), 22127 (for asthma), n_22040 // MET minutes per week


set more off
set maxvar 32000


use `pheno_data_1'

keep n_50_* ///  Height 
    n_23104_* /// BMI
    n_2050_* n_2060_* n_2070_* n_2080_* /// Depressive symptoms
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
    n_20107_* /// Father's illnesses
    n_20110_* /// Mother's illnesses
    n_1807_* /// Father's age at death
    n_3526_* /// Mother's age at death
    n_2946_* /// Father's age
    n_1845_* /// Mother's age
    /// n_3872* /// Age at First birth single child
    /// 20403_0_0  /// Amount of alcohol drunk on a typical drinking day
    /// 20405_0_0  /// Ever had known person concerned about, or recommend reduction of, alcohol consumption
    /// 20407_0_0  /// Frequency of failure to fulfil normal expectations due to drinking alcohol in last year
    /// 20408_0_0  /// Frequency of memory loss due to drinking alcohol in last year
    /// 20409_0_0  /// Frequency of feeling guilt or remorse after drinking alcohol in last year
    /// 20411_0_0  /// Ever been injured or injured someone else through drinking alcohol
    /// 20412_0_0  /// Frequency of needing morning drink of alcohol after heavy drinking session in last year
    /// 20413_0_0  /// Frequency of inability to cease drinking in last year
    /// 20414_0_0  /// Frequency of drinking alcohol
    /// 20416_0_0  /// Frequency of consuming six or more units of alcohol
    n_4559_* /// Family satisfaction
    n_4570_* /// Friendship satisfaction
    n_4537_* /// Work satisfaction
    n_4581_* /// Financial satisfaction
    n_2020_* /// loneliness
    n_3456_* n_2887_* n_6183_* /// Cigarettes per day
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
    n_6138_* /// EA
    n_20016_* /// IQ touchscreen
    n_20191_* /// IQ web-based
    n_2040_* /// Risk
    n_21000_* /// // Ethnicity
    n_21003_* /// age at assesment visit
    n_34_0_0 n_22200_0_0 /// birth year
    n_22010_* /// // Bad genotype n_22027_0_0
    n_eid // ID

gen IID=n_eid
merge 1:1 IID using `covar_data', nogen keep(match) 
merge 1:1 n_eid using `pheno_data_2', nogen keep(match)
merge 1:1 n_eid using `pheno_data_3', nogen keep(match)

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

forval i = 0/2 {
    gen AGE`i'=n_21003_`i'_0
    gen AGE`i'sq=AGE`i'*AGE`i'
    gen AGE`i'cb=AGE`i'*AGE`i'*AGE`i'
}

**********************************************************************************
******************************* GENERATE VARS ************************************
**********************************************************************************

**********************************************************
********************** ALZHEIMER'S  **********************
**********************************************************
* Father's ALZ status

foreach i in 0 1 2{
    gen ALZF_`i'=0 if n_20107_`i'_0!=.
    forval j = 0/3{
        replace ALZF_`i'=1 if n_20107_`i'_`j'==10 
    }
}


foreach i in 0 1 2{
    forval j = 0/3{
        replace ALZF_`i'=. if ALZF_`i'==0 & (n_20107_`i'_`j'==-11 | n_20107_`i'_`j'==-13)
    }
}

egen ALZF=rmax(ALZF_*)

foreach i in 0 1 2{
   gen F_age_`i'= n_2946_`i'_0 if n_2946_`i'_0>0
   replace F_age_`i' = n_1807_`i'_0 if F_age_`i'==. & n_1807_`i'_0>0
}

egen F_age=rmax(F_age_*) if ALZF==0
replace F_age=F_age_0 if ALZF_0==1
replace F_age=F_age_1 if ALZF_0==0 & ALZF_1==1
replace F_age=F_age_2 if ALZF_0==0 & ALZF_1==0 & ALZF_2==1

* Set to ALZF missing if father's ALZ status is 0 but father dies before age 60 or if father is younger than 60
replace ALZF=. if ALZF==0 & F_age<60

* Set to missing if father's age is not known
replace ALZF=. if F_age==.

**************************

* Mother's ALZ status
foreach i in 0 1 2{
    gen ALZM_`i'=0 if n_20110_`i'_0!=.
    forval j = 0/2{
        replace ALZM_`i'=1 if n_20110_`i'_`j'==10 
    }
}

foreach i in 0 1 2{
    forval j = 0/2{
        replace ALZM_`i'=. if ALZM_`i'==0 & (n_20110_`i'_`j'==-11 | n_20110_`i'_`j'==-13)
    }
}

egen ALZM=rmax(ALZM_*)

foreach i in 0 1 2{
   gen M_age_`i'= n_1845_`i'_0 if n_1845_`i'_0>0
   replace M_age_`i' = n_3526_`i'_0 if M_age_`i'==. & n_3526_`i'_0>0
}

egen M_age=rmax(M_age_*) if ALZM==0
replace M_age=M_age_0 if ALZM_0==1
replace M_age=M_age_1 if ALZM_0==0 & ALZM_1==1
replace M_age=M_age_2 if ALZM_0==0 & ALZM_1==0 & ALZM_2==1

* Set ALZM to missing if mother's ALZ status is 0 but mother dies before age 60 or if mother is younger than 60
replace ALZM=. if ALZM==0 & M_age<60

* Set to missing if father's age is not known
replace ALZM=. if M_age==.

**********************************************************

**********************************************************
************************ MENARCHE  ***********************
**********************************************************
forval i=0/3{
    replace n_2714_`i'_0=. if n_2714_`i'_0<0    
}
egen MENARCHE=rmean(n_2714_*)

**********************************************************



**********************************************************
************************* MYOPIA  ************************
**********************************************************

forval i=0/3{
    gen MYOPIA`i'=.
    replace MYOPIA`i'=0 if n_2207_`i'_0==0 | (n_6147_`i'_0!=. & n_6147_`i'_0>1)
    replace MYOPIA`i'=1 if n_6147_`i'_0==1
}
egen MYOPIA=rmax(MYOPIA*)

replace MYOPIA=0 if MYOPIA==. & n_20262_0_0==0
replace MYOPIA=1 if n_20262_0_0==1 | n_20262_0_0==2
  
**********************************************************

**********************************************************
******************* CANNABIS AMOUNT **********************
**********************************************************
ren n_20453_0_0 CANNABIS
replace CANNABIS=. if CANNABIS<0
**********************************************************


**********************************************************
******************** EVER CANNABIS  **********************
**********************************************************
gen EVERCANNABIS=0 if CANNABIS==0
replace EVERCANNABIS=1 if CANNABIS!=. & CANNABIS>0
**********************************************************


**********************************************************
*************************** AUDIT  ***********************
**********************************************************
foreach i in 403 405 407 408 409 411 412 413 414 416 {
    replace n_20`i'_0_0=. if n_20`i'_0_0<0
}


foreach i in 403 407 408 409 412 413 416{
    replace n_20`i'_0_0 = n_20`i'_0_0 - 1
}

foreach i in 405 411{
    replace n_20`i'_0_0 = 2*n_20`i'_0_0
}

egen AUDIT=rowtotal(n_20403_0_0 n_20405_0_0 n_20407_0_0 n_20408_0_0 n_20409_0_0 n_20411_0_0 n_20412_0_0 n_20413_0_0 n_20414_0_0 n_20416_0_0)

replace AUDIT=. if n_20414_0_0==.
replace AUDIT=. if n_20414_0_0==0 & (n_20405_0_0==. | n_20411_0_0==.)
replace AUDIT=. if n_20403_0_0 + n_20416_0_0>0 & (n_20407_0_0==. | n_20408_0_0==. | n_20409_0_0==. | n_20412_0_0==. | n_20413_0_0==. | n_20405_0_0==. | n_20411_0_0==.)
gen logAUDIT=log10(AUDIT)


**********************************************************
************************ EA ******************************
**********************************************************
forval i = 0/1 {
    forval j = 0/5 {
        g EA_`i'_`j' = 20 if n_6138_0_0 == 1
        replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 2
        replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 3
        replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 4
        replace EA_`i'_`j' = 19 if n_6138_`i'_`j' == 5
        replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
        replace EA_`i'_`j' = 7 if n_6138_`i'_`j' == -7
        replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
    }
}

*** Take max ***
egen EA = rmax(EA_*_*)
**********************************************************


**********************************************************
************************ IQ ******************************
**********************************************************

forval i = 0/2 {
    qui xi: reg n_20016_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict IQtouch`i',rstandard
}

qui xi: reg n_20191_0_0 SEX##c.AGE0 SEX##c.AGE0sq SEX##c.AGE0cb
predict IQweb, rstandard

egen IQ = rmean (IQtouch* IQweb)
**********************************************************



**********************************************************
************************ RISK ****************************
**********************************************************
forval i = 0/2 {
    replace n_2040_`i'_0=. if n_2040_`i'_0<0
    qui xi: reg n_2040_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict risk_`i', rstandard 
}

egen RISK=rmean(risk_*)

**********************************************************



**********************************************************
************************ HEIGHT **************************
**********************************************************
forval i = 0/2 {
    qui xi: reg n_50_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict height_`i', rstandard 
}

egen HEIGHT=rmean(height_*)
**********************************************************


**********************************************************
************************* BMI ****************************
**********************************************************
forval i = 0/1 {
    qui xi: reg n_23104_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict bmi_`i', rstandard 
}
egen BMI = rmean(bmi_*)
**********************************************************


**********************************************************
********************* NEUROTICISM ************************
**********************************************************
qui xi: reg n_20127_0_0 SEX##c.AGE0 SEX##c.AGE0sq SEX##c.AGE0cb
predict NEURO, rstandard
**********************************************************


**********************************************************
********************* MORNING PERSON *********************
**********************************************************
* scale 1-4, was reverse coded
forval i = 0/2 {
    gen morning_`i' =  5 - n_1180_`i'_0 if n_1180_`i'_0!=. & n_1180_`i'_0 > 0
    qui xi: reg morning_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict MORNING_`i', rstandard
}
egen MORNING = rmean(MORNING_*)
**********************************************************


**********************************************************
*********************** RELIGIOSITY **********************
**********************************************************
forval i = 0/2 {
    gen relig_`i' = 0 if (n_6160_`i'_0 != -3 & !missing(n_6160_`i'_0))
    forval j = 0/4 {
        replace relig_`i' = 1 if n_6160_`i'_`j' == 3
    }
    qui xi: reg relig_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict RELIG_`i', rstandard
}

egen RELIG = rmean(RELIG_*)
**********************************************************


**********************************************************
*************************** SWB **************************
**********************************************************
* scale 1-6, was reverse coded
forval i = 0/2 {
    gen swb_`i' = 7 - n_4526_`i'_0 if n_4526_`i'_0 > 0 
    qui xi: reg swb_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict SWB_`i', rstandard
}

egen SWB = rmean(SWB_*)
**********************************************************

*rerun
**********************************************************
*********************** DEPRESSION ***********************
**********************************************************

foreach i in 2050 2060 2070 2080{
    forval j=0/2{
        replace n_`i'_`j'_0=. if n_`i'_`j'_0<1
    }
}

forval i=0/2{
    gen depress_`i'=n_2050_`i'_0 + n_2060_`i'_0 + n_2070_`i'_0 + n_2080_`i'_0
    
    qui xi: reg depress_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict DEPRESS_`i', rstandard
}

egen DEPRESS = rmean(DEPRESS_*)
**********************************************************


**********************************************************
****************** FAMILY SATISFACTION *******************
**********************************************************
forval i=0/2 {
    gen famsat_`i'_0 = 7 - n_4559_`i'_0 if n_4559_`i'_0 > 0
    qui xi: reg famsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FAMSAT_`i', rstandard
}

egen FAM_SATISF = rmean(FAMSAT_*)
**********************************************************


**********************************************************
****************** FRIEND SATISFACTION *******************
**********************************************************
forval i=0/2 {
    gen friendsat_`i' = 7 - n_4570_`i'_0 if n_4570_`i'_0 > 0
    qui xi: reg friendsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FRIENDSAT_`i', rstandard
}

egen FRIEND_SATISF = rmean(FRIENDSAT_*)
**********************************************************


**********************************************************
***************** FINANCIAL SATISFACTION *****************
**********************************************************
forval i=0/2 {
    gen finsat_`i' = 7 - n_4581_`i'_0 if n_4581_`i'_0 > 0
    qui xi: reg finsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FINSAT_`i', rstandard
}

egen FIN_SATISF = rmean(FINSAT_*)
**********************************************************


**********************************************************
******************* WORK SATISFACTION ********************
**********************************************************
* 7 is unemployed
forval i=0/2 {
    gen worksat_`i' = 7 - n_4537_`i'_0 if (n_4537_`i'_0 > 0 & n_4537_`i'_0 != 7)
    qui xi: reg worksat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict WORKSAT_`i', rstandard
}

egen WORK_SATISF = rmean(WORKSAT_*)
**********************************************************


**********************************************************
********************** LONELINESS ************************
**********************************************************
forval i = 0/2 {
    replace n_2020_`i'_0=. if n_2020_`i'_0<0
    qui xi: reg n_2020_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict LONELY_`i', rstandard
}

egen LONELY = rmean(LONELY_*)
**********************************************************


**********************************************************
*********************** SELFHEALTH ***********************
**********************************************************
* scale 1-4, was reverse coded
forval i = 0/2 {
    gen health_`i' = 5 - n_2178_`i'_0 if n_2178_`i'_0 > 0 
    qui xi: reg health_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict HEALTH_`i', rstandard
}

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
    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
    (!missing(n_6152_2_0) & n_6152_2_0!=-3)) 


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

forval i = 0/3 {
    forval j = 0/2 {
        replace HAYFVR = 1 if n_6152_`j'_`i' == 9
    }
}
replace HAYFVR = 1 if n_6152_0_4 == 9


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
*gen ECZRHI = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
*    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
*    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
*    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
*    (!missing(n_6152_2_0) &  n_6152_2_0!=-3)) 


*forval i = 0/2 {
*    forval j = 0/3 {
*        replace ECZRHI = 1 if n_6152_`i'_`j' == 9
*    }
*}
    
*replace ECZRHI = 1 if n_6152_0_4 == 9

*forval i = 0/379 {
*    replace ECZRHI = 1 if (s_41202_0_`i' == "J301" | s_41202_0_`i' == "J302" | ///
*        s_41202_0_`i' == "J303" |  s_41202_0_`i' == "J304" | s_41202_0_`i' == "L20" | /// 
*        s_41202_0_`i' == "L208" |  s_41202_0_`i' == "J209")
*}

*forval i = 0/434 {
*    replace ECZRHI = 1 if (s_41204_0_`i' == "J301" | s_41204_0_`i' == "J302" | ///
*        s_41204_0_`i' == "J303" |  s_41204_0_`i' == "J304" | s_41204_0_`i' == "L20" | ///
*        s_41204_0_`i' == "L208" |  s_41204_0_`i' == "J209")
*}

*forval i = 0/28 {
*    replace ECZRHI = 1 if (n_20002_0_`i' == 1387 | n_20002_0_`i' == 1452)
*}

*forval i = 0/15 {
*    replace ECZRHI = 1 if (n_20002_1_`i' == 1387 | n_20002_1_`i' == 1452)
*}

*forval i = 0/16 {
*    replace ECZRHI = 1 if (n_20002_2_`i' == 1387 | n_20002_2_`i' == 1452)
*}
**********************************************************



**********************************************************
************** ASTHMA / HAYFEVER/ECZEMA ******************
**********************************************************
gen ASTECZRHI=.
replace ASTECZRHI=1 if (ASTHMA==1 | ECZEMA==1 | HAYFVR==1)
replace ASTECZRHI=0 if (ASTHMA==0 & ECZEMA==0 & HAYFVR==0)

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

    qui xi: reg AW_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict DPW_`i', rstandard
}

egen DPW = rmean(DPW_*)

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
keep n_eid SEX BATCH BYEAR IQ EA RISK HEIGHT BMI NEURO MORNING RELIG SWB HEALTH AFB N_EVERBORN_MEN  ///
    N_EVERBORN_WOMEN HAYFVR ECZEMA ASTECZRHI FAM_SATISF FRIEND_SATISF WORK_SATISF FIN_SATISF ///
    LONELY CPD CPDbins COPD ASTHMA MIGRAINE DEPRESS DPW EVERSMOKE IQ EA RISK CANNABIS EVERCANNABIS /// 
    logAUDIT MYOPIA MENARCHE ALZM ALZF F_age M_age PC*

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

* Save list of individuals in each partition
foreach partition in 1 2 3 {
    keep if partition == `partition'
    export delimited n_eid n_eid using "paritions/UKB_part`partition'_eid.txt", noq replace 
    clear
    use "tmp/pgs_repo.dta"
}
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
gen F_agesq=F_age*F_age
gen F_agecb=F_agesq*F_age
gen M_agesq=M_age*M_age
gen M_agecb=M_agesq*M_age

foreach partition in 1 2 3 {
    foreach var of varlist FAM_SATISF FRIEND_SATISF WORK_SATISF FIN_SATISF ///
        LONELY CPD CPDbins COPD ASTHMA MIGRAINE DEPRESS DPW EVERSMOKE ///
        HEIGHT BMI NEURO MORNING RELIG SWB HEALTH HAYFVR ECZEMA ECZRHI /// 
        ASTECZRHI EA IQ RISK CANNABIS EVERCANNABIS logAUDIT MYOPIA MENARCHE {
        
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


  
    qui xi:reg ALZF F_age F_agesq F_agecb i.BATCH PC1-PC40 if partition==`partition'
    predict resid, rstandard
    replace resid=. if partition!=`partition'
    export delimited FID IID resid using "input/UKB_ALZfather_part`partition'.pheno", noq delim(" ") replace

    drop resid

    qui xi:reg ALZM M_age M_agesq M_agecb i.BATCH PC1-PC40 if partition==`partition'
    predict resid, rstandard
    replace resid=. if partition!=`partition'
    export delimited FID IID resid using "input/UKB_ALZmother_part`partition'.pheno", noq delim(" ") replace

    drop resid


}


*** END ***
log close
