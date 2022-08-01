clear all

local WD="`1'"
local crosswalk="`2'"
local partition_data="`3'"
local pheno_data_1="`4'"
local pheno_data_2="`5'"
local pheno_data_3="`6'"
local pheno_data_4="`7'"
local covar_data="`8'"
local logfile="`9'"
local withdrawn="`10'"

cd `WD'
log using `logfile', replace
display "$S_DATE $S_TIME"

set more off
set maxvar 32000

use `pheno_data_1'

keep n_2139_* /// Age at first sex
    n_50_* ///  Height 
    n_23104_* /// BMI
    n_2050_* n_2060_* n_2070_* n_2080_* /// Depressive symptoms
    n_20127_0_0 /// Neuroticism score
    n_1180_* /// Morning person
    n_6160* /// Religious group activity
    n_4526*  /// Happiness / Subjective Wellbeing
    n_20160* /// Ever smoked
    n_2178* /// Self-rated Health
    n_2754* /// Age at First birth
    n_2734* /// Number of children (women)
    n_2405* /// Number of children (men)
    n_6152* /// Other diagnoses (Hayfever)
    n_6150* /// Other diagnoses (heart problems)
    n_1200_* /// Insomnia
    s_41202_* s_41204_* /// ICD10
    s_40001_* s_40002_* /// Death records - ICD10
    s_41203_* s_41205_* /// ICD9
    n_2453_* /// ever doctor diagnoed with cancer
    n_20001_* /// Self-reported cancer
    s_40006_* /// Type of cancer (ICD10)
    s_40013_* /// Type of cancer (ICD9)
    n_4079_* /// Diastolic blood pressure, automated reading
    n_94_* /// Diastolic blood pressure, manual reading
    n_4080_* /// Systolic blood pressure, automated reading
    n_93_* /// Systolic blood pressure, manual reading
    n_6153_* /// Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones	
    n_6177_* /// Medication for cholesterol, blood pressure or diabetes	
    n_20002_* /// Non-cancer illness
    ///n_3872* /// Age at first birth single child
    ///n_20403_0_0  /// Amount of alcohol drunk on a typical drinking day
    ///n_20405_0_0  /// Ever had known person concerned about, or recommend reduction of, alcohol consumption
    ///n_20407_0_0  /// Frequency of failure to fulfil normal expectations due to drinking alcohol in last year
    ///n_20408_0_0  /// Frequency of memory loss due to drinking alcohol in last year
    ///n_20409_0_0  /// Frequency of feeling guilt or remorse after drinking alcohol in last year
    ///n_20411_0_0  /// Ever been injured or injured someone else through drinking alcohol
    ///n_20412_0_0  /// Frequency of needing morning drink of alcohol after heavy drinking session in last year
    ///n_20413_0_0  /// Frequency of inability to cease drinking in last year
    ///n_20414_0_0  /// Frequency of drinking alcohol
    ///n_20416_0_0  /// Frequency of consuming six or more units of alcohol
    n_4559_* /// Family satisfaction
    n_4570_* /// Friendship satisfaction
    n_4537_* /// Work satisfaction
    n_4581_* /// Financial satisfaction
    n_2020_* /// Loneliness
    n_3456_* n_2887_* n_6183_* /// Cigarettes per day
    n_3436_* n_2867_* /// Age of smoking initiation
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
    n_20116_* /// smoking status
    n_1239_* /// current smoking
    n_1249_* /// past smoking
    n_2644_* /// light smoking
    n_6138_* /// EA
    n_845_* /// Age left schooling
    n_20016_* /// CP touchscreen
    n_20191_* /// CP web-based
    n_2040_* /// Risk
    n_21000_* /// // Ethnicity
    n_21003_* /// age at assesment visit
    n_34_0_0 n_22200_0_0 /// birth year
    n_22010_* /// // Bad genotype n_22027_0_0
    n_22001_* /// Genetic sex
    n_eid // ID

gen IID=n_eid
merge 1:1 IID using `covar_data', nogen keep(match) 
merge 1:1 n_eid using `pheno_data_2', nogen keep(match)
merge 1:1 n_eid using `pheno_data_3', nogen keep(match)
merge 1:1 n_eid using `pheno_data_4', nogen keep(match)

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

* Missing sex
drop if SEX == 0 
drop if missing(SEX)

* Sex mismatch
drop if InferredGender != SubmittedGender

* Missing batch
ren Batch BATCH
drop if missing(BATCH)

* Bad genotypes, initial release
drop if n_22010_0_0 == 1 
* Bad genotypes, second release
drop if hetmissingoutliers == 1
drop if missing(PC1-PC40)
drop if n_eid < 0 
drop if putativesexchromosomeaneuploidy == 1

* Drop non-Europeans
drop if PC1 > 0 

*** Keep White (value = 1) / British (value = 1001) / Irish (value = 1002) / Any other white background (value = 1003)  ***
gen ETHN =  n_21000_0_0
    replace ETHN = n_21000_1_0 if missing(ETHN)
    replace ETHN = n_21000_2_0 if missing(ETHN)
    drop if missing(ETHN)
    drop if (ETHN != 1 & ETHN != 1001 & ETHN != 1002 & ETHN != 1003)

*** Create BYEAR ***
gen BYEAR = (n_34_0_0 - 1900)/10
    replace BYEAR = (n_22200_0_0 - 1900)/10 if missing(BYEAR)
    drop if missing(BYEAR)

*** Create AGE ***
forval i = 0/3 {
    gen AGE`i'=n_21003_`i'_0
    gen AGE`i'sq=AGE`i'*AGE`i'
    gen AGE`i'cb=AGE`i'*AGE`i'*AGE`i'
}

**********************************************************************************
******************************* GENERATE VARS ************************************
**********************************************************************************

**********************************************************
******************** AGE AT FIRST BIRTH ******************
**********************************************************
* Minimum age of all times answered
forval i=0/2 {
    gen AFB1_`i' = n_2754_`i'_0 if n_2754_`i'_0 > 0
    gen AFB2_`i' = n_3872_`i'_0 if n_3872_`i'_0 > 0
}

egen AFB = rmin(AFB*)

**********************************************************

**********************************************************
******************** AGE AT FIRST SEX ********************
**********************************************************
* Minimum age of all times answered, <12 excluded as per Mills et al.
forval i=0/2 {
    gen AFS_`i' = n_2139_`i'_0 if n_2139_`i'_0 >= 12
}
egen AFS = rmin(AFS_*)
**********************************************************

**********************************************************
*************** AGE AT SMOKING INITIATION ****************
**********************************************************
* Minimum age of all times answered, <12 excluded as per Mills et al.
forval i=0/2 {
    gen ASI_current_`i' = n_3436_`i'_0 if n_3436_`i'_0 >= 0
    gen ASI_past_`i' = n_2867_`i'_0 if n_2867_`i'_0 >= 0
}
egen ASI = rmin(ASI_*)
**********************************************************

**********************************************************
******************** ANOREXIA NERVOSA ********************
**********************************************************
gen ANOREX=.

** CASES
* ICD10
forval i = 0/242 {
    replace ANOREX = 1 if s_41270_0_`i' == "F500" | s_41270_0_`i' == "F501"
}

* ICD9 
forval i = 0/46 {
    replace ANOREX = 1 if s_41271_0_`i' == "3071"
}

* Death register ICD10 - main
forval i = 0/2 { 
    replace ANOREX = 1 if s_40001_`i'_0 == "F500" |  s_40001_`i'_0 == "F501"
}

* Death register ICD10 - secondary
forval i = 0/14 {
    replace ANOREX = 1 if s_40002_0_`i' == "F500" | s_40002_0_`i' == "F501"
}
forval i = 0/9 {
    replace ANOREX = 1 if s_40002_1_`i' == "F500" | s_40002_1_`i' == "F501"
}
forval i = 0/7 {
    replace ANOREX = 1 if s_40002_2_`i' == "F500" | s_40002_2_`i' == "F501"
}

* Mental health problems
forval i=1/16{
    replace ANOREX = 1 if n_20544_0_`i' == 16
}

****************************

* CONTROLS assigned after BIPOLAR coding (using same set of controls as Bipolar and SCZ)
**********************************************************



**********************************************************
*********************** ALZHEIMER'S **********************
**********************************************************

**********************************************************


**********************************************************
************************** ASTHMA ************************
**********************************************************
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
replace AUDIT=logAUDIT

**********************************************************


**********************************************************
************************* BIPOLAR  ***********************
**********************************************************
gen BIPOLAR = 0

** CASES
* ICD10
forval i = 0/242 {
    forval j = 0/9 {    
        replace BIPOLAR = 1 if s_41270_0_`i' == "F31`j'"
    }
}

* ICD9 
forval i = 0/46 {
    foreach j in 2960 2961 2966 {
        replace BIPOLAR = 1 if s_41271_0_`i' == "`j'"
    }
}

* Self-report non-cancer illness
forval i = 0/28 {
    replace BIPOLAR = 1 if n_20002_0_`i' == 1291
}
forval i = 0/15 {
    replace BIPOLAR = 1 if n_20002_1_`i' == 1291
}
forval i = 0/33 {
    replace BIPOLAR = 1 if n_20002_2_`i' == 1291
}
forval i = 0/18 {
    replace BIPOLAR = 1 if n_20002_3_`i' == 1291
}


* Death register ICD10 - main
forval i = 0/2 {
    forval j = 0/9 {    
        replace BIPOLAR = 1 if s_40001_`i'_0 == "F31`j'"
    }
}

* Death register ICD10 - secondary
forval i = 0/14 {
    forval j = 0/9 {    
        replace BIPOLAR = 1 if s_40002_0_`i' == "F31`j'"
    }
}
forval i = 0/9 {
    forval j = 0/9 {    
        replace BIPOLAR = 1 if s_40002_1_`i' == "F31`j'"
    }
}
forval i = 0/7 {
    forval j = 0/9 {    
        replace BIPOLAR = 1 if s_40002_2_`i' == "F31`j'"
    }
}

* Bipolar disorder status
replace BIPOLAR = 1 if n_20122_0_0!=.

* Bipolar and major depression status
replace BIPOLAR = 1 if n_20126_0_0==1 | n_20126_0_0==2

* Source of report of F31 (bipolar affective disorder)
replace BIPOLAR = 1 if n_130893_0_0!=.

* Mental health problems
forval i=1/16{
    replace BIPOLAR = 1 if n_20544_0_`i' == 10
}

**********************************************

** CONTROLS: No report of mental health disorder in category 136 (Using only mental distress), did not report the use of any psychiatric medication at baseline (data-field 20003)
* Mental distress
*n_20499_0_0 == 0 //Ever sought or received professional help for mental distress
*n_20500_0_0 == 0 //Ever suffered mental distress preventing usual activities
*n_20544_0_0 == . //Mental health problems ever diagnosed by a professional

forval i=1/16{
    replace BIPOLAR = . if BIPOLAR == 0 & n_20544_0_`i' != .
}
replace BIPOLAR = . if BIPOLAR == 0 & (n_20499_0_0 != 0 | n_20500_0_0 != 0)

* Medications (https://eprints.gla.ac.uk/188130/2/188130Supp.pdf (eTable3))
forval i=0/47{
    * Mood stabilizers
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 {
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
    * Selective serotonin reuptake inhibitors
    foreach j in 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
    * Other antidepressants
    foreach j in 1141152732 1141152736 1141200564 1141201834 1141200570	1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
    * Traditional antipsychotics
    foreach j in 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
    * Second generation antipsychotics
    foreach j in 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
    * Sedatives & hypnotics
    foreach j in 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914	1140855920{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_0_`i'==`j'
    }
}
forval i=0/27{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_1_`i'==`j'
    }
}
forval i=0/29{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_2_`i'==`j'
    }
}
forval i=0/14{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace BIPOLAR=. if BIPOLAR==0 & n_20003_3_`i'==`j'
    }
}

* Depression
*n_20446_0_0 == 0 //Ever had prolonged feelings of sadness or depression
*n_20441_0_0 == 0 //Ever had prolonged loss of interest in normal activities
* Mania
*n_20502_0_0 == 0 //Ever had period extreme irritability
*n_20501_0_0 == 0 //Ever had period of mania / excitability
* Anxiety
*n_20421_0_0 == 0 //Ever felt worried, tense, or anxious for most of a month or longer
*n_20425_0_0 == 0 //Ever worried more than most people would in similar situation
* Unusual and psychotic experiences
*n_20468_0_0 == 0 //Ever believed in an un-real conspiracy against self
*n_20474_0_0 == 0 //Ever believed in un-real communications or signs
*n_20463_0_0 == 0 //Ever heard an un-real voice
*n_20471_0_0 == 0 //Ever seen an un-real vision
* Self-harm
*n_20480_0_0 == 0 //Ever self-harmed
* Addictions
*n_20401_0_0 == 0 //Ever addicted to any substance or behaviour

replace ANOREX = 0 if ANOREX ==. & BIPOLAR==0

**********************************************************


**********************************************************
****************** BLOOD PRESSURE ************************
**********************************************************
forval i = 0/3 {
    egen BPdia_`i'_auto = rmean(n_4079_`i'_0 n_4079_`i'_1)
    egen BPdia_`i'_manual = rmean(n_94_`i'_0 n_94_`i'_1) 
    egen BPdia_`i' = rmean(BPdia_`i'_auto BPdia_`i'_manual)
    replace BPdia_`i' = BPdia_`i' + 10  if n_6153_`i'_0 == 2 | n_6153_`i'_1 == 2 | n_6177_`i'_0 == 2 | n_6177_`i'_1 == 2
    replace BPdia_`i' = . if (n_6153_`i'_0 == -3 | n_6153_`i'_0 == -1 | n_6153_`i'_0 == .) & (n_6177_`i'_0 == -3 | n_6177_`i'_0 == -1 | n_6177_`i'_0 == .) 
    qui xi: reg BPdia_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict res_BPdia_`i', rstandard 
}
egen BPdia = rmean(res_BPdia_*)

forval i = 0/3 {
    egen BPsys_`i'_auto = rmean(n_4080_`i'_0 n_4080_`i'_1)
    egen BPsys_`i'_manual = rmean(n_93_`i'_0 n_93_`i'_1) 
    egen BPsys_`i' = rmean(BPsys_`i'_auto BPsys_`i'_manual)
    replace BPsys_`i' = BPsys_`i' + 15  if n_6153_`i'_0 == 2 | n_6153_`i'_1 == 2 | n_6177_`i'_0 == 2 | n_6177_`i'_1 == 2
    replace BPsys_`i' = . if (n_6153_`i'_0 == -3 | n_6153_`i'_0 == -1 | n_6153_`i'_0 == .) & (n_6177_`i'_0 == -3 | n_6177_`i'_0 == -1 | n_6177_`i'_0 == .) 
    qui xi: reg BPsys_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict res_BPsys_`i', rstandard 
}
egen BPsys = rmean(res_BPsys_*)

forval i = 0/3 {
    gen BPpulse_`i' = BPsys_`i' - BPdia_`i' 
    qui xi: reg BPpulse_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict res_BPpulse_`i', rstandard 
}
egen BPpulse = rmean(res_BPpulse_*)

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
******************** BREAST CANCER ***********************
**********************************************************
* ICD10 codes obtained from https://www.ambrygen.com/material/oncology/icd-10-code-reference-sheets/breast-cancer-icd-10-codes/630

* Diagnoses - ICD10
gen BRCA = 0

forval i = 0/242 {
    forval j = 500/509{
        replace BRCA = 1 if s_41270_0_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_41270_0_`i' == "`j'"
    }
}

**************************
* Diagnoses - ICD9 
forval i = 0/46 {
    forval j = 1740/1749 {
        replace BRCA = 1 if s_41271_0_`i' == "`j'"
    }
    foreach j in 1759 2330{
        replace BRCA = 1 if s_41271_0_`i' == "`j'"
    }
}

***************************
* Death register ICD10 - main
forval i = 0/2 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40001_`i'_0 == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40001_`i'_0 == "`j'"
    }
}


* Death register ICD10 - secondary
forval i = 0/14 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40002_0_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 0/9 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40002_1_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40002_1_`i' == "`j'"
    }
}
forval i = 0/7 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40002_2_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40002_2_`i' == "`j'"
    }
}

**************************

* Self-reported cancer
forval i = 0/5 {
    replace BRCA = 1 if n_20001_0_`i' == 1002
}
forval i = 0/3 {
    replace BRCA = 1 if n_20001_1_`i' == 1002
}
forval i = 0/4 {
    replace BRCA = 1 if n_20001_2_`i' == 1002
}
forval i = 0/3 {
    replace BRCA = 1 if n_20001_3_`i' == 1002
}

**************************
* Cancer register ICD10
forval i = 0/28 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40006_`i'_0 == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40006_`i'_0 == "`j'"
    }
}
forval i = 30/31 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40006_`i'_0 == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40006_`i'_0 == "`j'"
    }
}

**************************
* Cancer register ICD9
forval i = 0/13 {
    forval j = 1740/1749 {
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
    foreach j in 1759 2330{
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
}
foreach i in 15 16 17 18 19 22 25 26 27 29 {
    forval j = 1740/1749 {
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
    foreach j in 1759 2330{
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
}

replace BRCA = . if n_22001_0_0==1
**********************************************************


**********************************************************
*************************** CAD  *************************
**********************************************************

* HARDCAD CASES
gen HARDCAD=.

* ICD10
forval i = 0/242 {
    * Acute/subsequent myocardial infarction
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252" {
        replace HARDCAD = 1 if s_41270_0_`i' == "`j'"
    }
}

* ICD9 
forval i = 0/46 {
    foreach j in 4109 4119 4129 {
        replace HARDCAD = 1 if s_41271_0_`i' == "`j'"
    }
}


* Operations OPCS4
forval i = 0/123 {
    * Transluminal balloon angioplasty of coronary artery
    foreach j in "K491 K492 K493 K494 K498 K499" {
        replace HARDCAD = 1 if s_41272_0_`i' == "`j'"
    }
    * Other therapeutic transluminal operations on coronary artery
    foreach j in "K501 K502 K503 K504 K508 K509" {
        replace HARDCAD = 1 if s_41272_0_`i' == "`j'" 
    }
    * Off-pump coronary artery bypass grafting
    foreach j in "K401" "K402" "K403" "K404" "K408" "K409" "K411" "K412" "K413" "K414" "K418" "K419" "K421" "K422" "K423" "K424" "K428" "K429" "K431" "K432" "K433" "K434" "K438" "K439" "K441" "K442" "K448" "K449" "K451" "K452" "K453" "K454" "K455" "K456" "K458" "K459" "K461" "K462" "K463" "K464" "K465" "K468" "K469"{
        replace HARDCAD = 1 if s_41272_0_`i' == "`j'"
    }
    * Percutaneous transluminal balloon angioplasty and insertion of stent into coronary artery
    foreach j in "K751" "K752" "K753" "K754" "K758" "K759" {
        replace HARDCAD = 1 if s_41272_0_`i' == "`j'"
    }
}

* Operations OPCS3
forval i = 0/15 {
    foreach j in 3041 3042 3043 {
        replace HARDCAD = 1 if s_41273_0_`i' == "`j'"
    }
}

* Self-reported heart problems: heart attack
forval i = 0/2 {
    forval j = 0/3 {
        replace HARDCAD = 1 if n_6150_`i'_`j' == 1
    }
}

* Self-report non-cancer illness: heart attack
forval i = 0/28 {
    replace HARDCAD = 1 if n_20002_0_`i' == 1075
}
forval i = 0/15 {
    replace HARDCAD = 1 if n_20002_1_`i' == 1075
}
forval i = 0/33 {
    replace HARDCAD = 1 if n_20002_2_`i' == 1075
}
forval i = 0/18 {
    replace HARDCAD = 1 if n_20002_3_`i' == 1075
}

* Self report operations code
forval i = 0/31 {
    foreach j in 1070 1095 1523{
        replace HARDCAD = 1 if n_20004_0_`i' == `j'
    }
}
forval i = 0/14 {
    foreach j in 1070 1095 1523{
        replace HARDCAD = 1 if n_20004_1_`i' == `j'
    }
}
forval i = 0/17 {
    foreach j in 1070 1095 1523{
        replace HARDCAD = 1 if n_20004_2_`i' == `j'
    }
}
forval i = 0/9 {
    foreach j in 1070 1095 1523{
        replace HARDCAD = 1 if n_20004_3_`i' == `j'
    }
}


* Death register ICD10 - main
forval i = 0/2 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40001_`i'_0 == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40001_`i'_0 == "`j'"
    }
}

* Death register ICD10 - secondary
forval i = 0/14 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40002_0_`i' == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 0/9 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40002_1_`i' == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40002_1_`i' == "`j'"
    }
}
forval i = 0/7 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40002_2_`i' == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40002_2_`i' == "`j'"
    }
}


***********************************
* SOFTCAD - CASES
gen SOFTCAD = HARDCAD

* ICD10 - add angina and chronic ischaemic heart disease
forval i = 0/242 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259" {
        replace SOFTCAD = 1 if s_41270_0_`i' == "`j'"
    }
}


* ICD9  
forval i = 0/46 {
    foreach j in 4139 4140 4141 4148 4149{
        replace SOFTCAD = 1 if s_41271_0_`i' == "`j'"
    }
}


* Self-reported heart problems: angina
forval i = 0/2 {
    forval j = 0/3 {
        replace SOFTCAD = 1 if n_6150_`i'_`j' == 2
    }
}


* Self-report non-cancer illness: angina
forval i = 0/28 {
    replace SOFTCAD = 1 if n_20002_0_`i' == 1074
}
forval i = 0/15 {
    replace SOFTCAD = 1 if n_20002_1_`i' == 1074
}
forval i = 0/33 {
    replace SOFTCAD = 1 if n_20002_2_`i' == 1074
}
forval i = 0/18 {
    replace SOFTCAD = 1 if n_20002_3_`i' == 1074
}


* Death register ICD10 - main
forval i = 0/2 {
    * Angina / chronic ischaemic heart disease
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259" {
        replace SOFTCAD = 1 if s_40001_`i'_0 == "`j'" 
    }
}


* Death register ICD10 - secondary
forval i = 0/14 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259"{
        replace SOFTCAD = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 0/9 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259"{
        replace SOFTCAD = 1 if s_40002_1_`i' == "`j'"
    }
}
forval i = 0/7 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259"{
        replace SOFTCAD = 1 if s_40002_2_`i' == "`j'"
    }
}


***********************************

* CONTROLS
replace HARDCAD = 0 if SOFTCAD!=1 
replace SOFTCAD = 0 if SOFTCAD!=1

* Exclusions for aneurysm and atherosclerotic cardiovascular disease 
* ICD10
forval i = 0/242 {
    foreach j in "I250" "I253" "I254" {
        replace HARDCAD = . if s_41270_0_`i' == "`j'"
        replace SOFTCAD = . if s_41270_0_`i' == "`j'"
    }
}

* ICD10 death - main
forval i = 0/2 {
    foreach j in "I250" "I253" "I254" {
        replace HARDCAD = . if s_40001_`i'_0 == "`j'" 
        replace SOFTCAD = . if s_40001_`i'_0 == "`j'" 
    }
}

* ICD10 death - secondary
forval i = 0/14 {
    foreach j in "I250" "I253" "I254"{
        replace HARDCAD = . if s_40002_0_`i' == "`j'" 
        replace SOFTCAD = . if s_40002_0_`i' == "`j'" 
    }
}
forval i = 0/9 {
    foreach j in "I250" "I253" "I254"{
        replace HARDCAD = . if s_40002_1_`i' == "`j'" 
        replace SOFTCAD = . if s_40002_1_`i' == "`j'" 
    }
}
forval i = 0/7 {
    foreach j in "I250" "I253" "I254"{
        replace HARDCAD = . if s_40002_2_`i' == "`j'" 
        replace SOFTCAD = . if s_40002_2_`i' == "`j'" 
    }
}


* ICD9 main
forval i = 0/46 {
    foreach j in 4141 {
        replace HARDCAD = . if s_41271_0_`i' == "`j'"
        replace SOFTCAD = . if s_41271_0_`i' == "`j'"
    }
}

**********************************************************


**********************************************************
************************ CANNABIS  ***********************
**********************************************************
gen CANNABIS=.
replace CANNABIS=0 if n_20453_0_0==0
replace CANNABIS=1 if n_20453_0_0>0 & n_20453_0_0<=4
**********************************************************


**********************************************************
******************* CIGARETTES PER DAY *******************
**********************************************************
forval i = 0/2 {
    replace n_3456_`i'_0=. if n_3456_`i'_0 < 1 | n_3456_`i'_0 > 100
    replace n_2887_`i'_0=. if n_2887_`i'_0 < 1 | n_2887_`i'_0 > 100
    replace n_6183_`i'_0=. if n_6183_`i'_0 < 1 | n_6183_`i'_0 > 100
}

egen CPD = rmean(n_3456_* n_2887_* n_6183_*)
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


**********************************************************
**************** COGNITIVE PERFORMANCE *******************
**********************************************************
forval i = 0/2 {
    qui xi: reg n_20016_`i'_0 SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict CPtouch`i',rstandard
}

qui xi: reg n_20191_0_0 SEX##c.AGE0 SEX##c.AGE0sq SEX##c.AGE0cb
predict CPweb, rstandard

egen CP = rmean (CPtouch* CPweb)
**********************************************************


**********************************************************
****************** DEPRESSIVE SYMPTOMS *******************
**********************************************************
foreach i in 2050 2060 2070 2080{
    forval j=0/2{
        replace n_`i'_`j'_0=. if n_`i'_`j'_0<1
    }
}

forval i=0/2{
    gen depress_`i'=n_2050_`i'_0 + n_2060_`i'_0 + n_2070_`i'_0 + n_2080_`i'_0
    
    qui xi: reg depress_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict DEP_`i', rstandard
}

egen DEP = rmean(DEP_*)
**********************************************************


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
************************ EA ******************************
**********************************************************
forval i = 0/2 {
    replace n_845_`i'_0 = . if n_845_`i'_0 < 0
	forval j = 0/5 {
		g EA_`i'_`j' = 20 if n_6138_0_0 == 1
		replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 2
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 3
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 4
		*replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 5
        replace EA_`i'_`j' = n_845_`i'_0-5 if n_6138_`i'_`j' == 5
		replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
		replace EA_`i'_`j' = 7 if n_6138_`i'_`j' == -7
		replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
	}
}

*** Take max ***
egen EA = rmax(EA_*_*)
replace EA = . if EA < 7
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


**********************************************************
****************** FAMILY SATISFACTION *******************
**********************************************************
forval i=0/2 {
    gen famsat_`i'_0 = 7 - n_4559_`i'_0 if n_4559_`i'_0 > 0
    qui xi: reg famsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FAMSAT_`i', rstandard
}

egen FAMSAT = rmean(FAMSAT_*)
**********************************************************


**********************************************************
***************** FINANCIAL SATISFACTION *****************
**********************************************************
forval i=0/2 {
    gen finsat_`i' = 7 - n_4581_`i'_0 if n_4581_`i'_0 > 0
    qui xi: reg finsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FINSAT_`i', rstandard
}

egen FINSAT = rmean(FINSAT_*)
**********************************************************


**********************************************************
****************** FRIEND SATISFACTION *******************
**********************************************************
forval i=0/2 {
    gen friendsat_`i' = 7 - n_4570_`i'_0 if n_4570_`i'_0 > 0
    qui xi: reg friendsat_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict FRIENDSAT_`i', rstandard
}

egen FRIENDSAT = rmean(FRIENDSAT_*)
**********************************************************


**********************************************************
*********************** HAYFEVER *************************
**********************************************************
** Update when we get data field 22126

gen HAYFEVER = 0 if (!missing(s_41202_0_0) | !missing(s_41204_0_0) | ///
    !missing(n_20002_0_0) | !missing(n_20002_1_0) | !missing(n_20002_2_0) | ///
    (!missing(n_6152_0_0) & n_6152_0_0!=-3) | ///
    (!missing(n_6152_1_0) & n_6152_1_0!=-3) | ///
    (!missing(n_6152_2_0) & n_6152_2_0!=-3)) 


forval i = 0/379 {
    replace HAYFEVER = 1 if (s_41202_0_`i' == "J301" | s_41202_0_`i' == "J302" | ///
        s_41202_0_`i' == "J303" |  s_41202_0_`i' == "J304")
}

forval i = 0/434 {
    replace HAYFEVER = 1 if (s_41204_0_`i' == "J301" | s_41204_0_`i' == "J302" | ///
        s_41204_0_`i' == "J303" |  s_41204_0_`i' == "J304")
}

forval i = 0/28 {
    replace HAYFEVER = 1 if (n_20002_0_`i' == 1387 )
}

forval i = 0/15 {
    replace HAYFEVER = 1 if (n_20002_1_`i' == 1387 )
}

forval i = 0/16 {
    replace HAYFEVER = 1 if (n_20002_2_`i' == 1387 )
}

forval i = 0/3 {
    forval j = 0/2 {
        replace HAYFEVER = 1 if n_6152_`j'_`i' == 9
    }
}
replace HAYFEVER = 1 if n_6152_0_4 == 9
**********************************************************


**********************************************************
************** HAYFEVER / ASTHMA / ECZEMA ****************
**********************************************************
gen ASTECZRHI=.
replace ASTECZRHI=1 if (ASTHMA==1 | ECZEMA==1 | HAYFEVER==1)
replace ASTECZRHI=0 if (ASTHMA==0 & ECZEMA==0 & HAYFEVER==0)
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
************** INFLAMMATORY BOWEL DISEASE ****************
**********************************************************
gen IBD = 0 

* ICD10
forval i = 0/242 {
    forval j = 500/519 {
        replace IBD = 1 if s_41270_0_`i' == "K`j'"
    }
    foreach j in "K523" "K528" "K529" {
        replace IBD = 1 if s_41270_0_`i' == "`j'"
    }
}

* ICD9
forval i = 0/46 {
    foreach j in "5550" "5551" "5552" "5559" "5569" "5589" {
        replace IBD = 1 if s_41271_0_`i' == "`j'"
    }
}

* Death register ICD10 - main
forval i = 0/2 {
    forval j = 500/519 {
        replace IBD = 1 if s_40001_`i'_0 == "K`j'" 
    }
    foreach j in "K523" "K528" "K529" {
        replace IBD = 1 if s_40001_`i'_0 == "`j'"
    }    
}

* Death register ICD10 - secondary
forval i = 0/14 {
    forval j = 500/519 {
        replace IBD = 1 if s_40002_0_`i' == "K`j'" 
    }
    foreach j in "K523" "K528" "K529" {
        replace IBD = 1 if s_40002_0_`i' == "`j'"
    }    
}
forval i = 0/9 {
    forval j = 500/519 {
        replace IBD = 1 if s_40002_1_`i' == "K`j'" 
    }
    foreach j in "K523" "K528" "K529" {
        replace IBD = 1 if s_40002_1_`i' == "`j'"
    }    
}
forval i = 0/7 {
    forval j = 500/519 {
        replace IBD = 1 if s_40002_2_`i' == "K`j'" 
    }
    foreach j in "K523" "K528" "K529" {
        replace IBD = 1 if s_40002_2_`i' == "`j'"
    }    
}

**********************************************************


**********************************************************
************************ INSOMNIA ************************
**********************************************************

forval i = 0/2 {
    gen INSOMNIA_`i' = 0  if n_1200_`i'_0 != -3 & !missing(n_1200_`i'_0)
    replace INSOMNIA_`i' = 1  if n_1200_`i'_0 == 3
    qui xi: reg INSOMNIA_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict res_INSOMNIA_`i', rstandard
}
    
egen INSOMNIA=rmean(res_INSOMNIA_*)
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
************************ MENARCHE  ***********************
**********************************************************
forval i=0/3{
    replace n_2714_`i'_0=. if n_2714_`i'_0<0    
}
egen MENARCHE=rmean(n_2714_*)
**********************************************************


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
******************** NEARSIGHTEDNESS  ********************
**********************************************************
forval i=0/3{
    gen NEARSIGHTED`i'=.
    replace NEARSIGHTED`i'=0 if n_2207_`i'_0==0 | (n_6147_`i'_0!=. & n_6147_`i'_0>1)
    replace NEARSIGHTED`i'=1 if n_6147_`i'_0==1
}
egen NEARSIGHTED=rmax(NEARSIGHTED*)

replace NEARSIGHTED=0 if NEARSIGHTED==. & n_20262_0_0==0
replace NEARSIGHTED=1 if n_20262_0_0==1 | n_20262_0_0==2
**********************************************************


**********************************************************
********************* NEUROTICISM ************************
**********************************************************
qui xi: reg n_20127_0_0 SEX##c.AGE0 SEX##c.AGE0sq SEX##c.AGE0cb
predict NEURO, rstandard
**********************************************************


**********************************************************
******************** NUMBER EVER BORN ********************
**********************************************************
* max of number all reported children
gen NEBwomen_0 = n_2734_0_0 if n_2734_0_0 >= 0
gen NEBwomen_1 = n_2734_1_0 if n_2734_1_0 >= 0
gen NEBwomen_2 = n_2734_2_0 if n_2734_2_0 >= 0
egen NEBwomen = rmax(NEBwomen_*)

gen NEBmen_0 = n_2405_0_0 if n_2405_0_0 >= 0
gen NEBmen_1 = n_2405_1_0 if n_2405_1_0 >= 0
gen NEBmen_2 = n_2405_2_0 if n_2405_2_0 >= 0
egen NEBmen = rmax(NEBmen_*)
**********************************************************



**********************************************************
******************** PROSTATE CANCER *********************
**********************************************************
gen PRCA = 0

***************************
* Diagnoses - ICD10
forval i = 0/242 {
    replace PRCA = 1 if (s_41270_0_`i' == "C61" | s_41270_0_`i' == "D075")
}

**************************
* Diagnoses - ICD9 
forval i = 0/46 {
    replace PRCA = 1 if (s_41271_0_`i' == "1859" | s_41271_0_`i' == "2334")
}

***************************
* Death register ICD10 - main
forval i = 0/2 {
    replace PRCA = 1 if (s_40001_`i'_0 == "C61" | s_40001_`i'_0 == "D075") 
}


* Death register ICD10 - secondary
forval i = 0/14 {
    replace PRCA = 1 if (s_40002_0_`i' == "C61" | s_40002_0_`i' == "D075") 
}
forval i = 0/9 {
    replace PRCA = 1 if (s_40002_1_`i' == "C61" | s_40002_1_`i' == "D075") 
}
forval i = 0/7 {
    replace PRCA = 1 if (s_40002_2_`i' == "C61" | s_40002_2_`i' == "D075") 
}

***************************

* Self-reported cancer
forval i = 0/5 {
    replace PRCA = 1 if n_20001_0_`i' == 1044
}
forval i = 0/3 {
    replace PRCA = 1 if n_20001_1_`i' == 1044
}
forval i = 0/4 {
    replace PRCA = 1 if n_20001_2_`i' == 1044
}
forval i = 0/3 {
    replace PRCA = 1 if n_20001_3_`i' == 1044
}

**************************

* Cancer register ICD10
forval i = 0/28 {
    replace PRCA = 1 if (s_40006_`i'_0 == "C61" | s_40006_`i'_0 == "D075")
}

forval i = 30/31 {
    replace PRCA = 1 if (s_40006_`i'_0 == "C61" | s_40006_`i'_0 == "D075")
}

**************************
* Cancer register ICD9
forval i = 0/13 {
    replace PRCA = 1 if (s_40013_`i'_0 == "1859" | s_40013_`i'_0 == "2334")
}
foreach i in 15 16 17 18 19 22 25 26 27 29 {
    replace PRCA = 1 if (s_40013_`i'_0 == "1859" | s_40013_`i'_0 == "2334")
}

**************************

replace PRCA = . if n_22001_0_0==0
**********************************************************



**********************************************************
***************** RELIGIOUS ATTENDANCE *******************
**********************************************************
forval i = 0/2 {
    gen relig_`i' = 0 if (n_6160_`i'_0 != -3 & !missing(n_6160_`i'_0))
    forval j = 0/4 {
        replace relig_`i' = 1 if n_6160_`i'_`j' == 3
    }
    qui xi: reg relig_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict RELIGATT_`i', rstandard
}

egen RELIGATT = rmean(RELIGATT_*)
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
********************* SCHIZOPHRENIA **********************
**********************************************************
*gen SCZ = .
*replace SCZ = 1 if n_130875_0_0 != .
*forval i = 1/16 {
*    replace SCZ = 1 if n_20544_0_`i' == 2
*}

gen SCZ = .

** CASES: Schizophrenia / Schizoaffective disorders
* ICD10
forval i = 0/242 {
    forval j = 0/9 {    
        replace SCZ = 1 if s_41270_0_`i' == "F20`j'" |  s_41270_0_`i' == "F25`j'"
    }
}

* ICD9 
forval i = 0/46 {
    foreach j in 2953 2959 {
        replace SCZ = 1 if s_41271_0_`i' == "`j'"
    }
}

* Self-report non-cancer illness
forval i = 0/28 {
    replace SCZ = 1 if n_20002_0_`i' == 1289
}
forval i = 0/15 {
    replace SCZ = 1 if n_20002_1_`i' == 1289
}
forval i = 0/33 {
    replace SCZ = 1 if n_20002_2_`i' == 1289
}
forval i = 0/18 {
    replace SCZ = 1 if n_20002_3_`i' == 1289
}


* Death register ICD10 - main
forval i = 0/2 {
    forval j = 0/9 {    
        replace SCZ = 1 if s_40001_`i'_0 == "F20`j'" |  s_40001_`i'_0 == "F25`j'"
    }
}

* Death register ICD10 - secondary
forval i = 0/14 {
    forval j = 0/9 {    
        replace SCZ = 1 if s_40002_0_`i' == "F20`j'" |  s_40002_0_`i' == "F25`j'"
    }
}
forval i = 0/9 {
    forval j = 0/9 {    
        replace SCZ = 1 if s_40002_1_`i' == "F20`j'" |  s_40002_1_`i' == "F25`j'"
    }
}
forval i = 0/7 {
    forval j = 0/9 {    
        replace SCZ = 1 if s_40002_2_`i' == "F20`j'" |  s_40002_2_`i' == "F25`j'"
    }
}

* Source of report of F20/F25
replace SCZ = 1 if n_130875_0_0!=. | n_130885_0_0!=.

* Mental health problems
forval i=1/16{
    replace SCZ = 1 if n_20544_0_`i' == 2
}


** CONTROLS: Same as Bipolar controls
replace SCZ = 0 if SCZ==. & BIPOLAR==0
**********************************************************


**********************************************************
***************** SELF-RATED HEALTH **********************
**********************************************************
* scale 1-4, was reverse coded
forval i = 0/2 {
    gen health_`i' = 5 - n_2178_`i'_0 if n_2178_`i'_0 > 0 
    qui xi: reg health_`i' SEX##c.AGE`i' SEX##c.AGE`i'sq SEX##c.AGE`i'cb
    predict SELFHEALTH_`i', rstandard
}

egen SELFHEALTH = rmean(SELFHEALTH_*)
**********************************************************


**********************************************************
***************** SMOKING CESSATION **********************
**********************************************************
* scale 1-4, was reverse coded
gen SMCESS = .
forval i = 0/2 {
    gen SMCESS_`i' = 1 if n_20116_`i'_0==1
    replace SMCESS_`i' = 0 if n_20116_`i'_0==2
    replace SMCESS = SMCESS_`i' if SMCESS_`i'!=.
}
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


**********************************************************
********************* TYPE II DIABETES *******************
**********************************************************
gen T2D = 0 

* ICD10
forval i = 0/242 {
    forval j = 110/119 {
        replace T2D = 1 if s_41270_0_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D = . if T2D==0 & s_41270_0_`i' == "E`j'" 
    }
}

* ICD9
forval i = 0/46 {
    foreach j in "25000" "25010" "25020" "25090" {
        replace T2D = 1 if s_41271_0_`i' == "`j'" 
    }
    foreach j in "25009" "25019" "25029" "2503" "2504" "2505" "2506" "2507" "25099" {
        replace T2D = . if T2D==0 & s_41271_0_`i' == "`j'"
    }
}

* Death register ICD10 - main
forval i = 0/2 {
    forval j = 110/119 {
        replace T2D = 1 if s_40001_`i'_0 == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D = . if T2D==0 & s_40001_`i'_0 == "E`j'"
    }
}

* Death register ICD10 - secondary
forval i = 0/14 {
    forval j = 110/119 {
        replace T2D = 1 if s_40002_0_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D = . if T2D==0 & s_40002_0_`i' == "E`j'"
    }
}
forval i = 0/9 {
    forval j = 110/119 {
        replace T2D = 1 if s_40002_1_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D = . if T2D==0 & s_40002_1_`i' == "E`j'"
    }
}
forval i = 0/7 {
    forval j = 110/119 {
        replace T2D = 1 if s_40002_2_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D = . if T2D==0 & s_40002_2_`i' == "E`j'"
    }
}

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

egen WORKSAT = rmean(WORKSAT_*)
**********************************************************


**********************************************************************************
**********************************************************************************


*** SAVE FULL DATASET ***
keep n_eid SEX BATCH BYEAR AFB ANOREX ASTHMA AUDIT BIPOLAR BMI BPdia BPsys BPpulse ///
    BRCA HARDCAD SOFTCAD CANNABIS CP CPD COPD DEP DPW EA ECZEMA EVERSMOKE FAMSAT ///
    FINSAT FRIENDSAT HAYFEVER ASTECZRHI HEIGHT IBD LONELY MENARCHE MIGRAINE MORNING ///
    NEARSIGHTED NEBmen NEBwomen NEURO PRCA RELIGATT RISK SCZ SELFHEALTH SWB T2D WORKSAT PC*

save "tmp/pgi_repo.dta", replace


**********************************************************************************
**********************************************************************************

****************************************
**** Split into 3 equal partitions *****
****************************************
clear
import delimited `partition_data', delim("\t")
merge 1:1 n_eid using "tmp/pgi_repo.dta", keep(match using)

* Get maximum order in IDs_assignPartition_ordered.txt after merging 
qui sum partition_order

* For observations not in brain sample + relateds list, assign random value greater than maximum order to partition_order 
local maxorder =`r(max)'+1
replace partition_order=runiform()+`maxorder' if partition_order==.

* Assign into partitions based on partition order
sort partition_order
gen partition = ceil(3 * _n/_N)
save "tmp/pgi_repo.dta", replace

* Save list of individuals in each partition
foreach partition in 1 2 3 {
    keep if partition == `partition'
    export delimited n_eid n_eid using "partitions/UKB_part`partition'_eid.txt", noq delim(" ") replace 
    clear
    use "tmp/pgi_repo.dta"
}

********************************************************************************************************
********* Remove individuals that withdrew consent after the first release of the Repository ***********
***** Doing this here instead of in the beginning so that the partition assignment doesn't change ******
********************************************************************************************************
clear
import delimited `withdrawn'
ren v1 n_eid
merge 1:1 n_eid using "tmp/pgi_repo.dta", nogen keep(using)
save "tmp/pgi_repo.dta", replace


****************************************
********* Convert to GWAS IDs **********
****************************************
clear 
import delimited `crosswalk', delim(" ")
ren (v1 v2) (IID n_eid)
destring n_eid, replace ignore(" ")
merge 1:1 n_eid using "tmp/pgi_repo.dta", nogen
gen FID=IID


****************************************
********* RESIDUALIZE & EXPORT *********
****************************************

gen BYEAR2 = BYEAR*BYEAR
gen BYEAR3 = BYEAR2*BYEAR
gen SEXxBYEAR = SEX*BYEAR
gen SEXxBYEAR2 = SEX*BYEAR2
gen SEXxBYEAR3 = SEX*BYEAR3

foreach partition in 1 2 3 {
    foreach var of varlist BPdia BPsys BPpulse BIPOLAR HARDCAD SOFTCAD IBD T2D SCZ ANOREX ASTHMA //
        ASTECZRHI AUDIT BMI CANNABIS COPD CP CPD DEP DPW EA ECZEMA EVERSMOKE FAMSAT FINSAT FRIENDSAT ///
        HAYFEVER HEIGHT IBD LONELY MENARCHE MIGRAINE MORNING NEARSIGHTED NEURO RELIGATT RISK SELFHEALTH SWB WORKSAT {
        
        qui xi:reg `var' BYEAR* SEX* i.BATCH PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
        }


    * sex specific phenotypes:
    foreach var of varlist BRCA PRCA AFB NEBmen NEBwomen {
        
        qui xi:reg `var' BYEAR* i.BATCH PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }
    
}

*** END ***
log close
