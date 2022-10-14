clear all

local logfile="`1'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

use "tmp/pgi_repo.dta"


**********************************************************
*************** AGE AT SMOKING INITIATION ****************
**********************************************************
* Minimum age of all times answered
forval i=0/3 {
    gen ASI_current_`i' = n_3436_`i'_0 if n_3436_`i'_0 >= 0
    gen ASI_past_`i' = n_2867_`i'_0 if n_2867_`i'_0 >= 0
}
egen ASI = rmin(ASI_*)
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
************************ CANNABIS  ***********************
**********************************************************
gen CANNABIS=.
replace CANNABIS=0 if n_20453_0_0==0
replace CANNABIS=1 if n_20453_0_0>0 & n_20453_0_0<=4
**********************************************************


**********************************************************
******************* CIGARETTES PER DAY *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on max 51 individuals) only adds 34 individuals. Correlation is 0.9999, mean and SD identical to second decimal.
forval i = 0/3 {
    replace n_3456_`i'_0=. if n_3456_`i'_0 < 1 | n_3456_`i'_0 > 100
    replace n_2887_`i'_0=. if n_2887_`i'_0 < 1 | n_2887_`i'_0 > 100
    replace n_6183_`i'_0=. if n_6183_`i'_0 < 1 | n_6183_`i'_0 > 100
}

egen CPD = rmean(n_3456_* n_2887_* n_6183_*)

**********************************************************



**********************************************************
******************** DRINKS PER WEEK *********************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4870 individuals) only adds 6 individuals. Correlation is 0.9998, mean and SD identical to second decimal.
forval i = 0/3 {
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

forval i = 0/3 {
    * Set to sum of drinks if not less frequent than once or twice per week 
    egen AW_`i' = rowtotal(n_1568_`i'_0 n_1578_`i'_0 n_1588_`i'_0 n_1598_`i'_0 n_1608_`i'_0 n_5364_`i'_0) if  (n_1568_`i'_0!=. | n_1578_`i'_0!=. | n_1588_`i'_0!=. | n_1598_`i'_0!=. |  n_1608_`i'_0!=. |  n_5364_`i'_0!=.)

 
    * If less frequent than once or twice per week, replace with rescaled monthly value 
    egen AM_`i' = rowtotal(n_4407_`i'_0 n_4418_`i'_0 n_4429_`i'_0 n_4440_`i'_0 n_4451_`i'_0 n_4462_`i'_0) if (n_4407_`i'_0!=. | n_4418_`i'_0!=. | n_4429_`i'_0!=. | n_4440_`i'_0!=. | n_4451_`i'_0!=. | n_4462_`i'_0!=.)

    replace AW_`i'= AM_`i' / 4 if (n_1558_`i'_0 == 4 | n_1558_`i'_0 == 5 )
    
    * Set to 0 if never drinks 
    replace AW_`i'= 0 if n_1558_`i'_0 == 6

    qui xi: reg AW_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict DPW_`i', rstandard
}

egen DPW = rmean(DPW_*)
**********************************************************



**********************************************************
*********************** EVERSMOKE ************************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4869 individuals) adds no individuals. Correlation is 0.9999, mean and SD identical to third decimal
forval i = 0/3 {
    gen current_`i' = 1 if n_1239_`i'_0 == 1 
    replace current_`i' = 0 if n_1239_`i'_0 == 0 |  n_1239_`i'_0 == 2

    gen past_`i' = 1 if n_1249_`i'_0 == 1 | (n_1249_`i'_0 == 2 & n_2644_`i'_0 == 1) | (n_1249_`i'_0 == 3 & n_2644_`i'_0 == 1)
    replace past_`i' = 0 if n_1249_`i'_0 == 4 | (n_1249_`i'_0 == 2 & n_2644_`i'_0 == 0) | (n_1249_`i'_0 == 3 & n_2644_`i'_0 == 0)

    egen EVERSMOKE_`i'=rmax(current_`i' past_`i')
}

egen EVERSMOKE = rmax(EVERSMOKE_*)
**********************************************************


**********************************************************
***************** SMOKING CESSATION **********************
**********************************************************
gen SMCESS = .
forval i = 0/3 {
    gen SMCESS_`i' = 1 if n_20116_`i'_0==1
    replace SMCESS_`i' = 0 if n_20116_`i'_0==2
    replace SMCESS = SMCESS_`i' if SMCESS_`i'!=.
}
**********************************************************

**********************************************************************************
**********************************************************************************


*** SAVE FULL DATASET ***
keep n_eid FID IID Sex* Batch BYEAR ASI AUDIT CANNABIS CPD DPW EVERSMOKE SMCESS PC* partition

save "tmp/pgi_repo_substance.dta", replace


**********************************************************************************
**********************************************************************************

foreach partition in 1 2 3 {
    foreach var of varlist ASI AUDIT CANNABIS CPD DPW EVERSMOKE SMCESS {
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
        }    
}

*** END ***
log close
