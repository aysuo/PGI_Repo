clear all

local logfile="`1'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

use "tmp/pgi_repo.dta"
merge 1:1 n_eid using "tmp/gp_clinical.dta", nogen keep(master match)


**********************************************************
******************** PSYCH CONTROLS **********************
**********************************************************

** No report of mental health disorder in category 136 (Using only mental distress), did not report the use of any psychiatric medication at baseline (data-field 20003)
* Mental distress
*n_20499_0_0 == 0 //Ever sought or received professional help for mental distress
*n_20500_0_0 == 0 //Ever suffered mental distress preventing usual activities
*n_20544_0_0 == . //Mental health problems ever diagnosed by a professional

gen PsychCon = .
replace PsychCon = 0 if n_20499_0_0 == 0 & n_20500_0_0 == 0
forval i=1/16{
    replace PsychCon = . if n_20544_0_`i' != .
}


* Medications (https://eprints.gla.ac.uk/188130/2/188130Supp.pdf (eTable3))
forval i=0/47{
    * Mood stabilizers
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 {
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
    * Selective serotonin reuptake inhibitors
    foreach j in 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544{
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
    * Other antidepressants
    foreach j in 1141152732 1141152736 1141200564 1141201834 1141200570	1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812{
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
    * Traditional antipsychotics
    foreach j in 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312{
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
    * Second generation antipsychotics
    foreach j in 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320{
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
    * Sedatives & hypnotics
    foreach j in 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914	1140855920{
        replace PsychCon=. if n_20003_0_`i'==`j'
    }
}
forval i=0/27{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace PsychCon=. if n_20003_1_`i'==`j'
    }
}
forval i=0/29{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace PsychCon=. if n_20003_2_`i'==`j'
    }
}
forval i=0/14{
    foreach j in 1140867490 1140867504 1140867494 1140872198 1140872200 1141172838 1140872214 1140872064 2038459704 1140872072 1141167860 1141185460 1141162898 1140864452 1140872290 1140872302 1140867888	1140882236 1140879540 1140867876 1140921600 1141151946 1141180212 1141190158 1140867878 1140867884 1140879544 1141152732 1141152736 1141200564 1141201834 1141200570 1140916282 1140916288 1140879616 1140867658 1140867668 1140867662 1140867948 1140867934 1140867938 1140856186 1140867928 1140867850 1140910704 1140867852 1140867920 1140867922 1140879630 1140867712 1140867756 1140867758 1140879628 1140909806 1140867624 1141171824 1140879620 1140867690 1140867726 1140882310 1141146062 1140879556 1140867806 1140867812 1140879658 1140910358 1140863416 1140867168 1140867184 1140867092 1140867398 1140882098 1140867456 1140867156 1140856004 1140909800 1140867150 1140867152 1140867952 1140882100 1140867342 1140867406 1140867414 1140867084 1140867086 1140868120 1140867244 1140879750 1140867312 1141152848 1141152860 1140867444 1141177762 1140928916 1141167976 1141195974 1141202024 1141153490 1141184742 1140867420 1140882320 1140863152 1141157496 1140863244 1140863250 1140855856 1140863202 1140863210 1140863138 1140863144 1140928004 1141171404 1141171410 1140865016 1140864916 1140863182 1140863194 1140855896 1140863196 1140855900 1140855898 1140855902 1140855904 1140863104 1140863106 1140855914 1140855920{
        replace PsychCon=. if n_20003_3_`i'==`j'
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



**********************************************************
************************** ADHD  *************************
**********************************************************

** CASES
* ICD10
gen ADHD_ICD10 = .
forval i = 0/242 {
    forval j = 0/9 {    
        replace ADHD_ICD10 = 2 if s_41270_0_`i' == "F90`j'"
    }
}

* ICD9 
gen ADHD_ICD9 = .
forval i = 0/46 {
    forval j = 0/9 {  
        replace ADHD_ICD9 = 2 if s_41271_0_`i' == "314`j'"
    }
}

* Death register ICD10 - main
gen ADHD_ICD10Dm = .
forval i = 0/1 {
    forval j = 0/9 {    
        replace ADHD_ICD10Dm = 2 if s_40001_`i'_0 == "F90`j'"
    }
}

* Death register ICD10 - secondary
gen ADHD_ICD10Ds = .
forval i = 1/14 {
    forval j = 0/9 {    
        replace ADHD_ICD10Ds = 2 if s_40002_0_`i' == "F90`j'"
    }
}
forval i = 1/9 {
    forval j = 0/9 {    
        replace ADHD_ICD10Ds = 2 if s_40002_1_`i' == "F90`j'"
    }
}

* Mental health problems
gen ADHD_D = .
forval i=1/16{
    replace ADHD_D = 2 if n_20544_0_`i' == 18
}

replace ADHD_GP=. if ADHD_GP==0
egen ADHD=rmax(ADHD_ICD10 ADHD_ICD10Dm ADHD_ICD10Ds ADHD_ICD9 ADHD_D ADHD_GP PsychCon)
replace ADHD = . if ADHD==1
replace ADHD = 1 if ADHD==2

**********************************************************



**********************************************************
******************** ANOREXIA NERVOSA ********************
**********************************************************
** CASES
* ICD10
gen ANOREX_ICD10 = .
forval i = 0/242 {
    replace ANOREX_ICD10 = 2 if s_41270_0_`i' == "F500" | s_41270_0_`i' == "F501"
}

* ICD9 
gen ANOREX_ICD9 = .
forval i = 0/46 {
    replace ANOREX_ICD9 = 2 if s_41271_0_`i' == "3071"
}

* Death register ICD10 - main
gen ANOREX_ICD10Dm = .
forval i = 0/1 { 
    replace ANOREX_ICD10Dm = 2 if s_40001_`i'_0 == "F500" |  s_40001_`i'_0 == "F501"
}

* Death register ICD10 - secondary
gen ANOREX_ICD10Ds = .
forval i = 1/14 {
    replace ANOREX_ICD10Ds = 2 if s_40002_0_`i' == "F500" | s_40002_0_`i' == "F501"
}
forval i = 1/9 {
    replace ANOREX_ICD10Ds = 2 if s_40002_1_`i' == "F500" | s_40002_1_`i' == "F501"
}

* Mental health problems
gen ANOREX_SR = .
forval i=1/16{
    replace ANOREX_SR = 2 if n_20544_0_`i' == 16
}

replace ANOREX_GP=. if ANOREX_GP==0
egen ANOREX=rmax(ANOREX_ICD10 ANOREX_ICD10Dm ANOREX_ICD10Ds ANOREX_ICD9 ANOREX_SR ANOREX_GP PsychCon)
replace ANOREX = . if ANOREX==1
replace ANOREX = 1 if ANOREX==2

**********************************************************



**********************************************************
************************* BIPOLAR  ***********************
**********************************************************

** CASES
* ICD10
gen BIPOLAR_ICD10 = .
forval i = 0/242 {
    forval j = 0/9 {    
        replace BIPOLAR_ICD10 = 2 if s_41270_0_`i' == "F31`j'"
    }
}

* ICD9 
gen BIPOLAR_ICD9 = .
forval i = 0/46 {
    foreach j in 2960 2961 2966 {
        replace BIPOLAR_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}

* Self-report non-cancer illness
gen BIPOLAR_SR = .
forval i = 0/28 {
    replace BIPOLAR_SR = 2 if n_20002_0_`i' == 1291
}
forval i = 0/15 {
    replace BIPOLAR_SR = 2 if n_20002_1_`i' == 1291
}
forval i = 0/33 {
    replace BIPOLAR_SR = 2 if n_20002_2_`i' == 1291
}
forval i = 0/18 {
    replace BIPOLAR_SR = 2 if n_20002_3_`i' == 1291
}


* Death register ICD10 - main
gen BIPOLAR_ICD10Dm = .
forval i = 0/1 {
    forval j = 0/9 {    
        replace BIPOLAR_ICD10Dm = 2 if s_40001_`i'_0 == "F31`j'"
    }
}

* Death register ICD10 - secondary
gen BIPOLAR_ICD10Ds = .
forval i = 1/14 {
    forval j = 0/9 {    
        replace BIPOLAR_ICD10Ds = 2 if s_40002_0_`i' == "F31`j'"
    }
}
forval i = 1/9 {
    forval j = 0/9 {    
        replace BIPOLAR_ICD10Ds = 2 if s_40002_1_`i' == "F31`j'"
    }
}


* Bipolar disorder status
gen BIPOLAR_S1 = .
replace BIPOLAR_S1 = 2 if n_20122_0_0!=.

* Bipolar and major depression status
gen BIPOLAR_S2 = .
replace BIPOLAR_S2 = 2 if n_20126_0_0==1 | n_20126_0_0==2

* Source of report of F31 (bipolar affective disorder)
gen BIPOLAR_So = .
replace BIPOLAR_So = 2 if n_130893_0_0!=.

* Mental health problems
gen BIPOLAR_D = .
forval i=1/16{
    replace BIPOLAR_D = 2 if n_20544_0_`i' == 10
}

replace BIPOLAR_GP=. if BIPOLAR_GP==0
egen BIPOLAR=rmax(BIPOLAR_ICD10 BIPOLAR_ICD10Dm BIPOLAR_ICD10Ds BIPOLAR_ICD9 BIPOLAR_SR BIPOLAR_S1 BIPOLAR_S2 BIPOLAR_So BIPOLAR_D BIPOLAR_GP PsychCon)
replace BIPOLAR = . if BIPOLAR==1
replace BIPOLAR = 1 if BIPOLAR==2

**********************************************************


**********************************************************
****************** DEPRESSIVE SYMPTOMS *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4606 individuals) only adds 30 individuals. Correlation is 0.9997, mean and SD identical to third decimal.

foreach i in 2050 2060 2070 2080{
    forval j=0/3{
        replace n_`i'_`j'_0=. if n_`i'_`j'_0<1
    }
}

forval i=0/3{
    gen depress_`i'=n_2050_`i'_0 + n_2060_`i'_0 + n_2070_`i'_0 + n_2080_`i'_0
    
    qui xi: reg depress_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict DEP_`i', rstandard
}

egen DEP = rmean(DEP_*)
**********************************************************


**********************************************************
************************ INSOMNIA ************************
**********************************************************
forval i = 0/3 {
    gen INSOMNIA_`i' = 0  if n_1200_`i'_0 != -3 & !missing(n_1200_`i'_0)
    replace INSOMNIA_`i' = 1  if n_1200_`i'_0 == 3
    qui xi: reg INSOMNIA_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict res_INSOMNIA_`i', rstandard
}
    
egen INSOMNIA=rmean(res_INSOMNIA_*)
**********************************************************


**********************************************************
********************* SCHIZOPHRENIA **********************
**********************************************************

** CASES: Schizophrenia / Schizoaffective disorders
* ICD10
gen SCZ_ICD10 = .
forval i = 0/242 {
    forval j = 0/9 {    
        replace SCZ_ICD10 = 2 if s_41270_0_`i' == "F20`j'" |  s_41270_0_`i' == "F25`j'"
    }
}

* ICD9 
gen SCZ_ICD9 = .
forval i = 0/46 {
    foreach j in 2953 2959 {
        replace SCZ_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}

* Self-report non-cancer illness
gen SCZ_SR = .
forval i = 0/28 {
    replace SCZ_SR = 2 if n_20002_0_`i' == 1289
}
forval i = 0/15 {
    replace SCZ_SR = 2 if n_20002_1_`i' == 1289
}
forval i = 0/33 {
    replace SCZ_SR = 2 if n_20002_2_`i' == 1289
}
forval i = 0/18 {
    replace SCZ_SR = 2 if n_20002_3_`i' == 1289
}


* Death register ICD10 - main
gen SCZ_ICD10Dm = .
forval i = 0/1 {
    forval j = 0/9 {    
        replace SCZ_ICD10Dm = 2 if s_40001_`i'_0 == "F20`j'" |  s_40001_`i'_0 == "F25`j'"
    }
}

* Death register ICD10 - secondary
gen SCZ_ICD10Ds = .
forval i = 1/14 {
    forval j = 0/9 {    
        replace SCZ_ICD10Ds = 2 if s_40002_0_`i' == "F20`j'" |  s_40002_0_`i' == "F25`j'"
    }
}
forval i = 1/9 {
    forval j = 0/9 {    
        replace SCZ_ICD10Ds = 2 if s_40002_1_`i' == "F20`j'" |  s_40002_1_`i' == "F25`j'"
    }
}

* Source of report of F20/F25
gen SCZ_So = .
replace SCZ_So = 2 if n_130875_0_0!=. | n_130885_0_0!=.


* Mental health problems
gen SCZ_D = .
forval i=1/16{
    replace SCZ_D = 2 if n_20544_0_`i' == 2
}

replace SCZ_GP=. if SCZ_GP==0
egen SCZ=rmax(SCZ_ICD10 SCZ_ICD10Dm SCZ_ICD10Ds SCZ_ICD9 SCZ_SR SCZ_So SCZ_D SCZ_GP PsychCon)
replace SCZ = . if SCZ==1
replace SCZ = 1 if SCZ==2
**********************************************************

**********************************************************************************
**********************************************************************************


*** SAVE FULL DATASET ***
keep n_eid partition FID IID Sex Batch BYEAR ADHD ANOREX BIPOLAR DEP INSOMNIA SCZ PC*

save "tmp/pgi_repo_psych.dta", replace

sum ADHD ANOREX BIPOLAR DEP INSOMNIA SCZ

**********************************************************************************
**********************************************************************************

foreach partition in 1 2 3 {
    foreach var of varlist ADHD ANOREX BIPOLAR SCZ { /// DEP INSOMNIA {
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
        }
}

*** END ***
log close
