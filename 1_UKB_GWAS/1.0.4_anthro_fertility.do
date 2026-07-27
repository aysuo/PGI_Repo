clear all

local logfile="`1'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

use "tmp/pgi_repo.dta"


**********************************************************
******************** AGE AT FIRST BIRTH ******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on max 1606 individuals) adds no individuals. Correlation is 1.0000, mean and SD identical to third decimal.

* Minimum age of all times answered
forval i=0/3 {
    gen AFB1_`i' = n_2754_`i'_0 if n_2754_`i'_0 > 0
    gen AFB2_`i' = n_3872_`i'_0 if n_3872_`i'_0 > 0
}

egen AFB = rmin(AFB*)

**********************************************************

**********************************************************
******************** AGE AT FIRST SEX ********************
**********************************************************
* Minimum age of all times answered, <12 excluded as per Mills et al.
forval i=0/3 {
    gen AFS_`i' = n_2139_`i'_0 if n_2139_`i'_0 >= 12
}
egen AFS = rmin(AFS_*)
**********************************************************

**********************************************************
************************* BMI ****************************
**********************************************************
* Current GWAS excludes the third and fourth wave (as opposed to the code below). Not updating because adding the last waves (data on 47,201 and 4,772 individuals) adds 426 individuals. Correlation is 0.9976, mean and SD identical to second decimal.
forval i = 0/3 {
    qui xi: reg n_23104_`i'_0 Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict bmi_`i', rstandard 
}
egen BMI = rmean(bmi_*)
**********************************************************


**********************************************************
************************ HEIGHT **************************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the last wave (data on 4,881 individuals) adds no individuals. Correlation is 1.0000, mean and SD identical to third decimal
forval i = 0/3 {
    qui xi: reg n_50_`i'_0 Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict height_`i', rstandard 
}

egen HEIGHT=rmean(height_*)
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
*********** NUMBER EVER BORN, CHILDLESSNESS **************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the last wave (data on 2,506 individuals for women, 2,354 for men) adds no individuals for women and 1 for men. Correlation is 1.0000 for both men and women, mean and SD identical to third decimal

* max of number all reported children
forval i = 0/3 {
    gen NEBwomen_`i' = n_2734_`i'_0 if n_2734_`i'_0 >= 0 & AGE`i'>=45
    gen NEBmen_`i' = n_2405_`i'_0 if n_2405_`i'_0 >= 0 & AGE`i'>=55
}

egen NEBwomen = rmax(NEBwomen_*)
egen NEBmen = rmax(NEBmen_*)
gen CHILDLESSmen = 1 if  NEBmen==0
replace CHILDLESSmen = 0 if NEBmen>0 & NEBmen!=.
gen CHILDLESSwomen = 1 if NEBwomen==0
replace CHILDLESSwomen = 0 if NEBwomen>0 & NEBwomen!=.
**********************************************************


**********************************************************************************
**********************************************************************************


*** SAVE FULL DATASET ***
keep n_eid FID IID Sex* Batch BYEAR AFB AFS BMI HEIGHT MENARCHE CHILDLESSmen CHILDLESSwomen NEBmen NEBwomen PC* partition

save "tmp/pgi_repo_anthro_fertility.dta", replace


**********************************************************************************
**********************************************************************************
local anthro_fertility_sex_specific AFB MENARCHE NEBmen NEBwomen CHILDLESSmen CHILDLESSwomen
local anthro_fertility_pooled_sex AFS BMI HEIGHT
local anthro_fertility AFB AFS BMI HEIGHT MENARCHE NEBmen NEBwomen CHILDLESSmen CHILDLESSwomen
by partition, sort: summarize `anthro_fertility'


foreach partition in 1 2 3 {
    foreach var of varlist `anthro_fertility_pooled_sex' {
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }
    

    * sex specific phenotypes:
    foreach var of varlist `anthro_fertility_sex_specific' {
        qui xi:reg `var' BYEAR* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }
}


*** END ***
log close
