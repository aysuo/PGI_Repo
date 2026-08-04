clear all

local InputPhenoDir="`1'"
local OutputPhenoDir="`2'"
local logfile="`3'"

log using `logfile', replace
display "$S_DATE $S_TIME"

use "`InputPhenoDir'/pgi_repo_cog_person_wbeing.dta" 
merge 1:1 IID using "`InputPhenoDir'/pgi_repo_psych.dta", nogen 
merge 1:1 IID using "`InputPhenoDir'/pgi_repo_anthro_fertility.dta", nogen 
merge 1:1 IID using "`InputPhenoDir'/pgi_repo_health.dta", nogen 
merge 1:1 IID using "`InputPhenoDir'/pgi_repo_substance.dta", nogen 

drop FID IID
gen FID=n_eid
gen IID=n_eid
tabulate Batch, generate(batch)

save "`InputPhenoDir'/pgi_repo_prediction.dta", replace

foreach partition in 1 2 3 {
    use "`InputPhenoDir'/pgi_repo_prediction.dta", clear
    drop if partition!=`partition'
    export delimited FID IID PC1-PC20 batch* using "`OutputPhenoDir'/UKB`partition'/PC_BATCHdum.txt", noq delim(" ") replace
    export delimited FID IID PC1-PC20 batch* BYEAR* Sex* using "`OutputPhenoDir'/UKB`partition'/UKB`partition'.covar", noq delim(" ") replace
    foreach var of varlist ASI AUDIT CANNABIS CPD DPW EVERSMOKE SMCESS ALZ ALZnoproxy ALLERGYCAT ALLERGYPOLLEN ALLERGYDUST ASTHMA ASTECZRHI BL_LDL BL_HDL BL_CHOL BL_nonHDL BL_TRYG BPdia BPsys BPpulse COPD ECZEMA HARDCAD HAYFEVER IBD MIGRAINE NEARSIGHTED SOFTCAD T2D AFS BMI HEIGHT MENARCHE ADHD ANOREX ASD BIPOLAR SCZ DEP INSOMNIA CP EA FAMSAT FINSAT FRIENDSAT LONELY MORNING NEURO RELIGATT RISK SELFHEALTH SWB WORKSAT {       
        export delimited FID IID `var' using "`OutputPhenoDir'/UKB`partition'/`var'_noresid.pheno", noq delim(" ") replace
        qui xi:reg `var' BYEAR* Sex* 
        predict phenotype, rstandard
        export delimited FID IID phenotype using "`OutputPhenoDir'/UKB`partition'/`var'.pheno", noq delim(" ") replace
        drop phenotype
    }
}

foreach partition in 1 2 3 {
    use "`InputPhenoDir'/pgi_repo_prediction.dta", clear
    drop if partition!=`partition'
    foreach var of varlist AFB NEBmen NEBwomen BRCA PRCA CHILDLESSmen CHILDLESSwomen{    
        export delimited FID IID `var' using "`OutputPhenoDir'/UKB`partition'/`var'_noresid.pheno", noq delim(" ") replace
        qui xi:reg `var' BYEAR* 
        predict phenotype, rstandard
        export delimited FID IID phenotype using "`OutputPhenoDir'/UKB`partition'/`var'.pheno", noq delim(" ") replace
        drop phenotype
    }
}

log close

*** END ***

