clear all

local logfile="`1'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

use "tmp/pgi_repo.dta"


**********************************************************
**************** COGNITIVE PERFORMANCE *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4607 individuals) only adds 136 individuals. Correlation is 0.9991.
forval i = 0/2 {
    qui xi: reg n_20016_`i'_0 Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict CPtouch`i',rstandard
}

qui xi: reg n_20191_0_0 Sex##c.AGE0 Sex##c.AGE0sq Sex##c.AGE0cb
predict CPweb, rstandard

egen CP = rmean (CPtouch* CPweb)
**********************************************************


**********************************************************
************************ EA ******************************
**********************************************************
egen ageLeftSchool = rmax(n_845_*_0) 
replace ageLeftSchool = . if ageLeftSchool < 0

forval i = 0/3 {
	forval j = 0/5 {
		gen EA_`i'_`j' = 20 if n_6138_`i'_`j' == 1
		replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 2
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 3
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 4
        replace EA_`i'_`j' = ageLeftSchool-5 if n_6138_`i'_`j' == 5
		replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
		replace EA_`i'_`j' = 7 if n_6138_`i'_`j' == -7
		replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
	}
}
egen EA = rmax(EA_*_*)
replace EA = . if EA < 7
**********************************************************


**********************************************************
****************** FAMILY SATISFACTION *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4835 individuals) only adds 21 individuals. Correlation is 0.9982.
forval i=0/3 {
    gen famsat_`i'_0 = 7 - n_4559_`i'_0 if n_4559_`i'_0 > 0
    qui xi: reg famsat_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict FAMSAT_`i', rstandard
}

egen FAMSAT = rmean(FAMSAT_*)
**********************************************************


**********************************************************
***************** FINANCIAL SATISFACTION *****************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4864 individuals) only adds 20 individuals. Correlation is 0.9984.
forval i=0/3 {
    gen finsat_`i' = 7 - n_4581_`i'_0 if n_4581_`i'_0 > 0
    qui xi: reg finsat_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict FINSAT_`i', rstandard
}

egen FINSAT = rmean(FINSAT_*)
**********************************************************


**********************************************************
****************** FRIEND SATISFACTION *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4820 individuals) only adds 37 individuals. Correlation is 0.9982.
forval i=0/3 {
    gen friendsat_`i' = 7 - n_4570_`i'_0 if n_4570_`i'_0 > 0
    qui xi: reg friendsat_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict FRIENDSAT_`i', rstandard
}

egen FRIENDSAT = rmean(FRIENDSAT_*)
**********************************************************


**********************************************************
********************** LONELINESS ************************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave  (data on 4796 individuals) only adds 3 individuals. Correlation is 0.9996.
forval i = 0/3 {
    replace n_2020_`i'_0=. if n_2020_`i'_0<0
    qui xi: reg n_2020_`i'_0 Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict LONELY_`i', rstandard
}

egen LONELY = rmean(LONELY_*)
**********************************************************


**********************************************************
********************* MORNING PERSON *********************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4454 individuals) only adds 83 individuals. Correlation is 0.9998.
* scale 1-4, was reverse coded
forval i = 0/3 {
    gen morning_`i' =  5 - n_1180_`i'_0 if n_1180_`i'_0!=. & n_1180_`i'_0 > 0
    qui xi: reg morning_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict MORNING_`i', rstandard
}
egen MORNING = rmean(MORNING_*)
**********************************************************


**********************************************************
********************* NEUROTICISM ************************
**********************************************************
qui xi: reg n_20127_0_0 Sex##c.AGE0 Sex##c.AGE0sq Sex##c.AGE0cb
predict NEURO, rstandard
**********************************************************


**********************************************************
***************** RELIGIOUS ATTENDANCE *******************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4860 individuals) adds no individuals. Correlation is 0.9999.
forval i = 0/3 {
    gen relig_`i' = 0 if (n_6160_`i'_0 != -3 & !missing(n_6160_`i'_0))
    forval j = 0/4 {
        replace relig_`i' = 1 if n_6160_`i'_`j' == 3
    }
    qui xi: reg relig_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict RELIGATT_`i', rstandard
}

egen RELIGATT = rmean(RELIGATT_*)
**********************************************************


**********************************************************
************************ RISK ****************************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4666 individuals) only adds 8 individuals. Correlation is 0.9996.
forval i = 0/3 {
    replace n_2040_`i'_0=. if n_2040_`i'_0<0
    qui xi: reg n_2040_`i'_0 Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict risk_`i', rstandard 
}

egen RISK=rmean(risk_*)
**********************************************************


**********************************************************
***************** SELF-RATED HEALTH **********************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4866 individuals) adds no individuals. Correlation is 0.9996.
* scale 1-4, was reverse coded
forval i = 0/3 {
    gen health_`i' = 5 - n_2178_`i'_0 if n_2178_`i'_0 > 0 
    qui xi: reg health_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict SELFHEALTH_`i', rstandard
}

egen SELFHEALTH = rmean(SELFHEALTH_*)
**********************************************************


**********************************************************
*************************** SWB **************************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 4862 individuals) only adds 27 individuals. Correlation is 0.9984.
* scale 1-6, was reverse coded
forval i = 0/3 {
    gen swb_`i' = 7 - n_4526_`i'_0 if n_4526_`i'_0 > 0 
    qui xi: reg swb_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict SWB_`i', rstandard
}

egen SWB = rmean(SWB_*)
**********************************************************


**********************************************************
******************* WORK SATISFACTION ********************
**********************************************************
* Current GWAS excludes the fourth wave (as opposed to the code below). Not updating because adding the fourth wave (data on 2713 individuals) only adds 325 individuals. Correlation is 0.9982.
* 7 is unemployed
forval i=0/3 {
    gen worksat_`i' = 7 - n_4537_`i'_0 if (n_4537_`i'_0 > 0 & n_4537_`i'_0 != 7)
    qui xi: reg worksat_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict WORKSAT_`i', rstandard
}

egen WORKSAT = rmean(WORKSAT_*)
**********************************************************



*** SAVE FULL DATASET ***
keep n_eid FID IID Sex* Batch BYEAR partition CP EA FAMSAT FINSAT FRIENDSAT LONELY MORNING NEURO RELIGATT RISK SELFHEALTH SWB WORKSAT PC*

save "tmp/pgi_repo_cog_person_wbeing.dta", replace


**********************************************************************************
**********************************************************************************

foreach partition in 1 2 3 {
    foreach var of varlist EA { /// CP EA FAMSAT FINSAT FRIENDSAT LONELY MORNING NEURO RELIGATT RISK SELFHEALTH SWB WORKSAT {
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
        }
}

*** END ***
log close
