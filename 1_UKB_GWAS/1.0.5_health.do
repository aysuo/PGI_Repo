clear all

local logfile="`1'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

use "tmp/pgi_repo.dta"
merge 1:1 n_eid using "tmp/gp_clinical.dta", nogen keep(master match)



**********************************************************
********************** ALZHEIMER'S  **********************
**********************************************************
* Father's ALZ status
forval i = 0/3{
    gen ALZF_`i' = 0 if n_20107_`i'_0 != . & n_20107_`i'_0 != -11 & n_20107_`i'_0 != -13 
    forval j = 0/3{
        replace ALZF_`i'=1 if n_20107_`i'_`j'==10 
    }
}

forval i = 0/3{
    gen F_age_`i'= n_2946_`i'_0 if n_2946_`i'_0!=. & n_2946_`i'_0>0 
    replace F_age_`i' = n_1807_`i'_0 if  n_1807_`i'_0!=. & n_1807_`i'_0>0
}


* Mother's ALZ status
forval i = 0/3{
    gen ALZM_`i'=0 if n_20110_`i'_0 != . & n_20110_`i'_0!= -11 & n_20110_`i'_0!= -13
    forval j = 0/3{
        replace ALZM_`i'=1 if n_20110_`i'_`j'==10 
    }
}

forval i = 0/3{
    gen M_age_`i'= n_1845_`i'_0 if n_1845_`i'_0!=. & n_1845_`i'_0>0
    replace M_age_`i' = n_3526_`i'_0 if n_3526_`i'_0!=. & n_3526_`i'_0>0
}


foreach i in M F {
    egen ALZ`i'=rmax(ALZ`i'_*)

    * Latest age with non-missing ALZ status
    gen `i'_age = `i'_age_3 if ALZ`i'_3 != .
    replace `i'_age = `i'_age_2 if ALZ`i'_3 == . & ALZ`i'_2 != .
    replace `i'_age = `i'_age_1 if ALZ`i'_3 == . & ALZ`i'_2 == . & ALZ`i'_1 != .
    replace `i'_age = `i'_age_0 if ALZ`i'_3 == . & ALZ`i'_2 == . & ALZ`i'_1 == . & ALZ`i'_0 != .

    replace ALZ`i' = (100-`i'_age) / 100 if ALZ`i'==0
    replace ALZ`i' = 0.32 if ALZ`i' > 0.32 & ALZ`i' < 1
    replace ALZ`i' = 0 if ALZ`i' < 0
}

gen ALZ_P = ALZF + ALZM 

**************************

* In-patient - ICD10
gen ALZ_ICD10 = .
forval i = 0/242 {
    foreach j in "F000" "F001" "F002" "F009" "G300" "G301" "G308" "G309"{
        replace ALZ_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}


* Death register ICD10 - main
gen ALZ_ICD10Dm = .
forval i = 0/1 {
    foreach j in "F000" "F001" "F002" "F009" "G300" "G301" "G308" "G309" {
        replace ALZ_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}


* Death register ICD10 - secondary
gen ALZ_ICD10Ds = .
forval i = 1/14 {
    foreach j in "F000" "F001" "F002" "F009" "G300" "G301" "G308" "G309" {
        replace ALZ_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "F000" "F001" "F002" "F009" "G300" "G301" "G308" "G309" {
        replace ALZ_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}


* In-patient - ICD9
gen ALZ_ICD9 = .
forval i = 0/46 {
    foreach j in "3310" {
        replace ALZ_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}


replace ALZ_GP=. if ALZ_GP==0 
egen ALZi=rmax(ALZ_ICD10 ALZ_ICD10Dm ALZ_ICD10Ds ALZ_ICD9 ALZ_GP)

gen ALZnoproxy = 0
replace ALZnoproxy = . if ALZi==1
replace ALZnoproxy = 1 if ALZi==2

replace ALZi = . if ALZi==1
egen ALZ=rmax(ALZi ALZ_P)


**********************************************************

**********************************************************


**********************************************************
*******************  ALLERGY - CAT ***********************
**********************************************************
ren ALLERGYCAT_GP ALLERGYCAT
replace ALLERGYCAT = . if ALLERGYCAT==1
replace ALLERGYCAT = 1 if ALLERGYCAT==2
**********************************************************

**********************************************************
*******************  ALLERGY - DUST ***********************
**********************************************************
* Self report non-cancer illness
gen ALLERGYDUST_SR = 0
forval i = 0/28 {
    replace ALLERGYDUST_SR = 2 if (n_20002_0_`i' == 1668 )
}

forval i = 0/15 {
    replace ALLERGYDUST_SR = 2 if (n_20002_1_`i' == 1668 )
}

forval i = 0/33 {
    replace ALLERGYDUST_SR = 2 if (n_20002_2_`i' == 1668 )
}

forval i = 0/18 {
    replace ALLERGYDUST_SR = 2 if (n_20002_3_`i' == 1668 )
}
**************************
egen ALLERGYDUST=rmax(ALLERGYDUST_SR ALLERGYDUST_GP)
replace ALLERGYDUST = . if ALLERGYDUST==1
replace ALLERGYDUST = 1 if ALLERGYDUST==2
**********************************************************


**********************************************************
*******************  ALLERGY - POLLEN ********************
**********************************************************
gen ALLERGYPOLLEN_ICD10 = 0
* In-patient - ICD10
forval i = 0/242 {
    foreach j in "J301" {
        replace ALLERGYPOLLEN_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
* Exclude from controls: Allergic rhinitis - unspecified, predominantly allergic asthma
forval i = 0/242 {
    foreach j in "J304" "J450" {
        replace ALLERGYPOLLEN_ICD10 = 1 if ALLERGYPOLLEN_ICD10 == 0 & s_41270_0_`i' == "`j'"
    }
}

**************************

* Death register ICD10 - main
gen ALLERGYPOLLEN_ICD10Dm=0
forval i = 0/1 {
    foreach j in "J301" {
        replace ALLERGYPOLLEN_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}
forval i = 0/1 {
    foreach j in "J304" "J450" {
        replace ALLERGYPOLLEN_ICD10Dm = 1 if ALLERGYPOLLEN_ICD10Dm == 0 & s_40001_`i'_0 == "`j'"
    }
}


**************************

* Death register ICD10 - secondary
gen ALLERGYPOLLEN_ICD10Ds=0
forval i = 1/14 {
    foreach j in "J301" {
        replace ALLERGYPOLLEN_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J301" {
        replace ALLERGYPOLLEN_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}
forval i = 1/14 {
    foreach j in "J304" "J450" {
        replace ALLERGYPOLLEN_ICD10Ds = 1 if ALLERGYPOLLEN_ICD10Ds == 0 & s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J304" "J450" {
        replace ALLERGYPOLLEN_ICD10Ds = 1 if ALLERGYPOLLEN_ICD10Ds == 0 & s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9
gen ALLERGYPOLLEN_ICD9 = 0
forval i = 0/46 {
    foreach j in "4770" {
        replace ALLERGYPOLLEN_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}
forval i = 0/46 {
    foreach j in "4779" {
        replace ALLERGYPOLLEN_ICD9 = 1 if ALLERGYPOLLEN_ICD9 == 0 & s_41271_0_`i' == "`j'"
    }
}

**************************
egen ALLERGYPOLLEN=rmax(ALLERGYPOLLEN_ICD10 ALLERGYPOLLEN_ICD10Dm ALLERGYPOLLEN_ICD10Ds ALLERGYPOLLEN_ICD9 ALLERGYPOLLEN_GP)
replace ALLERGYPOLLEN = . if ALLERGYPOLLEN==1
replace ALLERGYPOLLEN = 1 if ALLERGYPOLLEN==2
**********************************************************


**********************************************************
************************** ASTHMA ************************
**********************************************************
gen ASTHMA_ICD10 = 0
* In-patient - ICD10
forval i = 0/242 {
    foreach j in "J450" "J451" "J458" "J459" "J46" {
        replace ASTHMA_ICD10 = 1 if s_41270_0_`i' == "`j'"
    }
}

**************************

* Death register ICD10 - main
gen ASTHMA_ICD10Dm=0
forval i = 0/1 {
    foreach j in "J450" "J451" "J458" "J459" "J46" {
        replace ASTHMA_ICD10Dm = 1 if s_40001_`i'_0 == "`j'"
    }
}

**************************

* Death register ICD10 - secondary
gen ASTHMA_ICD10Ds=0
forval i = 1/14 {
    foreach j in "J450" "J451" "J458" "J459" "J46" {
        replace ASTHMA_ICD10Ds = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J450" "J451" "J458" "J459" "J46" {
        replace ASTHMA_ICD10Ds = 1 if s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9
gen ASTHMA_ICD9 = 0
forval i = 0/46 {
    foreach j in "49300" "49309" "49310" "49319" "49390" "49399" {
        replace ASTHMA_ICD9 = 1 if s_41271_0_`i' == "`j'"
    }
}

**************************

* Self report non-cancer illness
gen ASTHMA_SR = 0
forval i = 0/28 {
    replace ASTHMA_SR = 1 if (n_20002_0_`i' == 1111 )
}

forval i = 0/15 {
    replace ASTHMA_SR = 1 if (n_20002_1_`i' == 1111 )
}

forval i = 0/33 {
    replace ASTHMA_SR = 1 if (n_20002_2_`i' == 1111 )
}

forval i = 0/18 {
    replace ASTHMA_SR = 1 if (n_20002_3_`i' == 1111 )
}

**************************

* Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor	
gen ASTHMA_D = 0
forval i = 0/3 {
    forval j = 0/2 {
        replace ASTHMA_D = 1 if n_6152_`j'_`i' == 8
    }
}
replace ASTHMA_D = 1 if n_6152_3_0 == 8 | n_6152_3_1 == 8 | n_6152_3_2 == 8

**************************
* "Has a doctor ever told you that you have had any of the conditions below?" "asthma"
gen ASTHMA_D2 = 0
replace ASTHMA_D2 = 1 if n_22127_0_0 == 1

**************************

replace ASTHMA_GP=1 if ASTHMA_GP==2
egen ASTHMA=rmax(ASTHMA_ICD10 ASTHMA_ICD10Dm ASTHMA_ICD10Ds ASTHMA_ICD9 ASTHMA_SR ASTHMA_D ASTHMA_D2 ASTHMA_GP)

**********************************************************


**********************************************************
******************** BLOOD LIPIDS ************************
**********************************************************
* Statin adjustment factors
*                   ID      Multiplier	Offset	P
*LDL                30780	1.46155	    0	    6.87E-218
*Total cholesterol  30690	1.33588	    0	    4.66E-213


gen n_statins_0 = 0
forval i=0/47{
    foreach j in 1140861958 1140888594 1140888648 1141146234 1141192410 1140861922 1141146138{
        replace n_statins_0 = n_statins_0 + 1 if n_20003_0_`i'==`j'
    }
}
gen n_statins_1 = 0
forval i=0/27{
    foreach j in 1140861958 1140888594 1140888648 1141146234 1141192410 1140861922 1141146138{
        replace n_statins_1 = n_statins_1 + 1 if n_20003_1_`i'==`j'
    }
}

forval i=0/1{
    gen BL_LDL_`i' = n_30780_`i'_0
    gen BL_HDL_`i' = n_30760_`i'_0
    gen BL_CHOL_`i' = n_30690_`i'_0
    gen BL_TRYG_`i' = n_30870_`i'_0

    * Statin adjustment
    replace BL_LDL_`i' = n_30780_`i'_0 * 1.46155 if n_statins_`i' > 0
    replace BL_CHOL_`i' = n_30690_`i'_0 * 1.33588 if n_statins_`i' > 0

    gen BL_nonHDL_`i' = BL_CHOL_`i'-BL_HDL_`i'

    * Age indicator
    gen AGEind_`i' = AGE`i'
    replace AGEind_`i' = 50 if AGE`i'<50
    replace AGEind_`i' = 78 if AGE`i'>78

    * 5-year age intervals
    gen AGEbin_`i' = floor(AGE`i'/5)

    * Icosatiles of sampling time
    xtile stime_`i'= ts_3166_`i'_0, nq(20)

    * Fasting time
    gen ftime_`i' = n_74_`i'_0
    replace ftime_`i' = 18 if n_74_`i'_0 > 18
    replace ftime_`i' = 1 if n_74_`i'_0 == 0

    * Icosatiles of estimated sample dilution factor
    xtile dilution_`i' = n_30897_`i'_0, nq(20)

    * Assesment center
    ren n_54_`i'_0 acenter_`i'

    * We don't have assessment date, approximate from byear,age at assessment and assessment month
    gen amonth_`i' = 100*(n_34_0_0 + AGE`i') + n_55_`i'_0

    * Aliquot
    gen al_LDL_`i' = n_30782_`i'_0
    gen al_HDL_`i' = n_30762_`i'_0
    gen al_TRYG_`i' = n_30872_`i'_0
    gen al_CHOL_`i' = n_30692_`i'_0

    * Paper also includes assay time for each measure but we didn't get those variables (30781, 30761, 30871, 30691)
}

* Truncate assessment month at the ends (too few individuals)
replace amonth_0 = 200500 if amonth_0 < 200605
replace amonth_0 = 201008 if amonth_0 >= 201008 & amonth_0 <= 201010
replace amonth_1 = 201100 if amonth_1 < 201109


foreach j in LDL HDL CHOL TRYG{
    forval i=0/1{
        gen logBL_`j'_`i'=log10(BL_`j'_`i')
        qui xi: reg logBL_`j'_`i' i.AGEind_`i' Sex  Sex#i.AGEbin_`i' i.stime_`i' i.ftime_`i' i.dilution_`i' i.acenter_`i' i.amonth_`i' i.al_`j'_`i' 
        predict res_BL_`j'_`i', rstandard 
    }
    egen BL_`j' = rmean(res_BL_`j'_*)
}

forval i=0/1{
    gen logBL_nonHDL_`i'=log10(BL_nonHDL_`i')
    qui xi: reg logBL_nonHDL_`i' i.AGEind_`i' Sex  Sex#i.AGEbin_`i' i.stime_`i' i.ftime_`i' i.dilution_`i' i.acenter_`i' i.amonth_`i' i.al_HDL_`i' i.al_CHOL_`i'
    predict res_BL_nonHDL_`i', rstandard 
}
egen BL_nonHDL = rmean(res_BL_nonHDL_*)

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
    qui xi: reg BPdia_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict res_BPdia_`i', rstandard 
}
egen BPdia = rmean(res_BPdia_*)

forval i = 0/3 {
    egen BPsys_`i'_auto = rmean(n_4080_`i'_0 n_4080_`i'_1)
    egen BPsys_`i'_manual = rmean(n_93_`i'_0 n_93_`i'_1) 
    egen BPsys_`i' = rmean(BPsys_`i'_auto BPsys_`i'_manual)
    replace BPsys_`i' = BPsys_`i' + 15  if n_6153_`i'_0 == 2 | n_6153_`i'_1 == 2 | n_6177_`i'_0 == 2 | n_6177_`i'_1 == 2
    replace BPsys_`i' = . if (n_6153_`i'_0 == -3 | n_6153_`i'_0 == -1 | n_6153_`i'_0 == .) & (n_6177_`i'_0 == -3 | n_6177_`i'_0 == -1 | n_6177_`i'_0 == .) 
    qui xi: reg BPsys_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict res_BPsys_`i', rstandard 
}
egen BPsys = rmean(res_BPsys_*)

forval i = 0/3 {
    gen BPpulse_`i' = BPsys_`i' - BPdia_`i' 
    qui xi: reg BPpulse_`i' Sex##c.AGE`i' Sex##c.AGE`i'sq Sex##c.AGE`i'cb
    predict res_BPpulse_`i', rstandard 
}
egen BPpulse = rmean(res_BPpulse_*)

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
forval i = 0/1 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40001_`i'_0 == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40001_`i'_0 == "`j'"
    }
}


* Death register ICD10 - secondary
forval i = 1/14 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40002_0_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40002_1_`i' == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40002_1_`i' == "`j'"
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
forval i = 0/21 {
    forval j = 500/509 {
        replace BRCA = 1 if s_40006_`i'_0 == "C`j'"
    }
    foreach j in "D050" "D051" "D057" "D059" {
        replace BRCA = 1 if s_40006_`i'_0 == "`j'"
    }
}

**************************
* Cancer register ICD9
forval i = 0/8 {
    forval j = 1740/1749 {
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
    foreach j in 1759 2330{
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
}
forval i = 10/12 {
    forval j = 1740/1749 {
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
    foreach j in 1759 2330{
        replace BRCA = 1 if s_40013_`i'_0 == "`j'"
    }
}
foreach i in 14 {
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
forval i = 0/1 {
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
forval i = 1/14 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40002_0_`i' == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    * MI
    foreach j in "I210" "I211" "I212" "I213" "I214" "I219" "I220" "I221" "I228" "I229" "I230" "I231" "I232" "I233" "I234" "I235" "I236" "I238" "I240" "I241" "I248" "I249" "I252"{
        replace HARDCAD = 1 if s_40002_1_`i' == "`j'"
    }
    * Complications related to coronary bypass/angioplasty/graft
    foreach j in "T822" "Z951" "Z955" {
        replace HARDCAD = 1 if s_40002_1_`i' == "`j'"
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
    foreach j in 4139 4140 4141 4148 4149 {
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
forval i = 0/1 {
    * Angina / chronic ischaemic heart disease
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259" {
        replace SOFTCAD = 1 if s_40001_`i'_0 == "`j'" 
    }
}


* Death register ICD10 - secondary
forval i = 1/14 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259"{
        replace SOFTCAD = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "I200" "I201" "I208" "I209" "I250" "I251" "I253" "I254" "I255" "I256" "I258" "I259"{
        replace SOFTCAD = 1 if s_40002_1_`i' == "`j'"
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
forval i = 0/1 {
    foreach j in "I250" "I253" "I254" {
        replace HARDCAD = . if s_40001_`i'_0 == "`j'" 
        replace SOFTCAD = . if s_40001_`i'_0 == "`j'" 
    }
}

* ICD10 death - secondary
forval i = 1/14 {
    foreach j in "I250" "I253" "I254"{
        replace HARDCAD = . if s_40002_0_`i' == "`j'" 
        replace SOFTCAD = . if s_40002_0_`i' == "`j'" 
    }
}
forval i = 1/9 {
    foreach j in "I250" "I253" "I254"{
        replace HARDCAD = . if s_40002_1_`i' == "`j'" 
        replace SOFTCAD = . if s_40002_1_`i' == "`j'" 
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
************************** COPD **************************
**********************************************************
* ICD10 codes received from here, interstitial/compensatory emphysema excluded
* https://www.cigna.com/static/docs/starplus/icd10-copd.pdf

gen COPD_ICD10 = 0
* In-patient - ICD10
forval i = 0/242 {
    foreach j in "J440" "J441" "J448" "J449" "J410" "J411" "J418" "J42" "J430" "J431" "J432" "J438" "J439"{
        replace COPD_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
forval i = 0/242 {
    foreach j in "J40"{
        replace COPD_ICD10 = 1 if COPD_ICD10 == 0 & s_41270_0_`i' == "`j'"
    }
}
**************************

* Death register ICD10 - main
gen COPD_ICD10Dm=0
forval i = 0/1 {
    foreach j in "J440" "J441" "J448" "J449" "J410" "J411" "J418" "J42" "J430" "J431" "J432" "J438" "J439" {
        replace COPD_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}
forval i = 0/1 {
    foreach j in "J40" {
        replace COPD_ICD10Dm = 1 if COPD_ICD10Dm == 0 & s_40001_`i'_0 == "`j'"
    }
}

**************************

* Death register ICD10 - secondary
gen COPD_ICD10Ds=0
forval i = 1/14 {
    foreach j in "J440" "J441" "J448" "J449" "J410" "J411" "J418" "J42" "J430" "J431" "J432" "J438" "J439" {
        replace COPD_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J440" "J441" "J448" "J449" "J410" "J411" "J418" "J42" "J430" "J431" "J432" "J438" "J439" {
        replace COPD_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}
forval i = 1/14 {
    foreach j in "J40" {
        replace COPD_ICD10Ds = 1 if COPD_ICD10Ds == 0 & s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J40" {
        replace COPD_ICD10Ds = 1 if COPD_ICD10Ds == 0 & s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9
gen COPD_ICD9 = 0
forval i = 0/46 {
    foreach j in "4910" "4911" "4912" "4918" "4919" "4929" "4969"{
        replace COPD_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}
forval i = 0/46 {
    foreach j in "4909"{
        replace COPD_ICD9 = 1 if COPD_ICD9 == 0 & s_41271_0_`i' == "`j'"
    }
}


**************************
* Self report non-cancer illness
gen COPD_SR = 0
forval i = 0/28 {
    replace COPD_SR = 2 if (n_20002_0_`i' == 1112 | n_20002_0_`i' == 1113)
}

forval i = 0/15 {
    replace COPD_SR = 2 if (n_20002_1_`i' == 1112 | n_20002_1_`i' == 1113)
}

forval i = 0/33 {
    replace COPD_SR = 2 if (n_20002_2_`i' == 1112 | n_20002_2_`i' == 1113)
}

forval i = 0/18 {
    replace COPD_SR = 2 if (n_20002_3_`i' == 1112 | n_20002_3_`i' == 1113)
}

**************************
* Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor	
gen COPD_D = 0
forval i = 0/3 {
    forval j = 0/2 {
        replace COPD_D = 2 if n_6152_`j'_`i' == 6
    }
}
replace COPD_D = 2 if n_6152_3_0 == 6 | n_6152_3_1 == 6 | n_6152_3_2 == 6


**************************
* "Has a doctor ever told you that you have had any of the conditions below?" "emphysema" (22128),"chronic bronchitis" (22129), "COPD" (22130)
gen COPD_D2 = 0
replace COPD_D2 = 2 if n_22128_0_0 == 1 | n_22129_0_0 == 1 |  n_22130_0_0 == 1

**************************
egen COPD=rmax(COPD_ICD10 COPD_ICD10Dm COPD_ICD10Ds COPD_ICD9 COPD_SR COPD_D COPD_D2 COPD_GP)
replace COPD = . if COPD==1
replace COPD = 1 if COPD==2
**********************************************************




**********************************************************
************************ ECZEMA **************************
**********************************************************
gen ECZEMA_ICD10 = 0
* In-patient - ICD10
forval i = 0/242 {
    foreach j in "L20" "L208" "L209"{
        replace ECZEMA_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
* Exclude from controls: Exfoliative dermatitis, lichen simplex chronicus, prurigo, pruritus infective dermatitis, pityriasis alba, unspecified dermatitits as these could be forms of atopic dermatitis, conditions following it, or symptoms of it.
forval i = 0/242 {
    foreach j in "L26" "L28" "L280" "L281" "L282" "L29" "L290" "L291" "L292" "L293" "L298" "L299" "L303" "L305" "L309" {
        replace ECZEMA_ICD10 = 1 if ECZEMA_ICD10 == 0 & s_41270_0_`i' == "`j'"
    }
}
**************************

* Death register ICD10 - main
gen ECZEMA_ICD10Dm=0
forval i = 0/1 {
    foreach j in "L20" "L208" "L209" {
        replace ECZEMA_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}
forval i = 0/1 {
    foreach j in "L26" "L28" "L280" "L281" "L282" "L29" "L290" "L291" "L292" "L293" "L298" "L299" "L303" "L305" "L309" {
        replace ECZEMA_ICD10Dm = 1 if ECZEMA_ICD10Dm == 0 & s_40001_`i'_0 == "`j'"
    }
}

**************************

* Death register ICD10 - secondary
gen ECZEMA_ICD10Ds=0
forval i = 1/14 {
    foreach j in "L20" "L208" "L209" {
        replace ECZEMA_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "L20" "L208" "L209" {
        replace ECZEMA_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}
forval i = 1/14 {
    foreach j in "L26" "L28" "L280" "L281" "L282" "L29" "L290" "L291" "L292" "L293" "L298" "L299" "L303" "L305" "L309" {
        replace ECZEMA_ICD10Ds = 1 if ECZEMA_ICD10Ds == 0 & s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "L26" "L28" "L280" "L281" "L282" "L29" "L290" "L291" "L292" "L293" "L298" "L299" "L303" "L305" "L309" {
        replace ECZEMA_ICD10Ds = 1 if ECZEMA_ICD10Ds == 0 & s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9
gen ECZEMA_ICD9 = 0
forval i = 0/46 {
    foreach j in "691" "6918" "69180" {
        replace ECZEMA_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}
forval i = 0/46 {
    foreach j in "690" "6909" "697" "6970" "6971" "6978" "6979" "698" "6980" "6981" "6982" "6983" "6988" "6989"{
        replace ECZEMA_ICD9 = 1 if ECZEMA_ICD9 == 0 & s_41271_0_`i' == "`j'"
    }
}

**************************

* Self report non-cancer illness
* Too broad (eczema/dermatitis), so not including in cases but removing from controls

gen ECZEMA_SR = 0
forval i = 0/28 {
    replace ECZEMA_SR = 1 if (n_20002_0_`i' == 1452 )
}

forval i = 0/15 {
    replace ECZEMA_SR = 1 if (n_20002_1_`i' == 1452 )
}

forval i = 0/33 {
    replace ECZEMA_SR = 1 if (n_20002_2_`i' == 1452 )
}

forval i = 0/18 {
    replace ECZEMA_SR = 1 if (n_20002_3_`i' == 1452 )
}

**************************
* Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor	
* Too broad (hayfever / allergic rhinitis / eczema), so not including in cases but removing from controls

gen ECZEMA_D = 0
forval i = 0/3 {
    forval j = 0/2 {
        replace ECZEMA_D = 1 if n_6152_`j'_`i' == 9
    }
}
replace ECZEMA_D = 1 if n_6152_3_0 == 9 | n_6152_3_1 == 9 | n_6152_3_2 == 9

**************************
**************************
egen ECZEMA=rmax(ECZEMA_ICD10 ECZEMA_ICD10Dm ECZEMA_ICD10Ds ECZEMA_ICD9 ECZEMA_SR ECZEMA_D ECZEMA_GP)
replace ECZEMA = . if ECZEMA==1
replace ECZEMA = 1 if ECZEMA==2
**********************************************************



**********************************************************
*********************** HAYFEVER *************************
**********************************************************

gen HAYFEVER_ICD10 = 0
* In-patient - ICD10
forval i = 0/242 {
    foreach j in "J301" "J302" "J303" "J304"{
        replace HAYFEVER_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
* Exclude from controls: Vasomotor and allergic rhinitis
forval i = 0/242 {
    foreach j in "J30" {
        replace HAYFEVER_ICD10 = 1 if HAYFEVER_ICD10 == 0 & s_41270_0_`i' == "`j'"
    }
}
**************************
* Death register ICD10 - main
gen HAYFEVER_ICD10Dm=0
forval i = 0/1 {
    foreach j in "J301" "J302" "J303" "J304" {
        replace HAYFEVER_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}
forval i = 0/1 {
    foreach j in "J30" {
        replace HAYFEVER_ICD10Dm = 1 if HAYFEVER_ICD10Dm == 0 & s_40001_`i'_0 == "`j'"
    }
}

**************************
* Death register ICD10 - secondary
gen HAYFEVER_ICD10Ds=0
forval i = 1/14 {
    foreach j in "J301" "J302" "J303" "J304" {
        replace HAYFEVER_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J301" "J302" "J303" "J304" {
        replace HAYFEVER_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}
forval i = 1/14 {
    foreach j in "J30" {
        replace HAYFEVER_ICD10Ds = 1 if HAYFEVER_ICD10Ds == 0 & s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "J30" {
        replace HAYFEVER_ICD10Ds = 1 if HAYFEVER_ICD10Ds == 0 & s_40002_1_`i' == "`j'"
    }
}

**************************
* In-patient - ICD9
gen HAYFEVER_ICD9 = 0
forval i = 0/46 {
    foreach j in "477" "4770" "4778" "4779" {
        replace HAYFEVER_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}

**************************
* Self report non-cancer illness

gen HAYFEVER_SR = 0
forval i = 0/28 {
    replace HAYFEVER_SR = 2 if (n_20002_0_`i' == 1387 )
}

forval i = 0/15 {
    replace HAYFEVER_SR = 2 if (n_20002_1_`i' == 1387 )
}

forval i = 0/33 {
    replace HAYFEVER_SR = 2 if (n_20002_2_`i' == 1387 )
}

forval i = 0/18 {
    replace HAYFEVER_SR = 2 if (n_20002_3_`i' == 1387 )
}

**************************
* Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor	
* Too broad (hayfever / allergic rhinitis / eczema), so not including in cases but removing from controls

gen HAYFEVER_D = 0
forval i = 0/3 {
    forval j = 0/2 {
        replace HAYFEVER_D = 1 if n_6152_`j'_`i' == 9
    }
}
replace HAYFEVER_D = 1 if n_6152_3_0 == 9 | n_6152_3_1 == 9 | n_6152_3_2 == 9

**************************
* "Has a doctor ever told you that you have had any of the conditions below?" "hayfever or allergic rhinitis"
gen HAYFEVER_D2 = 0
replace HAYFEVER_D2 = 2 if n_22126_0_0 == 1

**************************
**************************

egen HAYFEVER=rmax(HAYFEVER_ICD10 HAYFEVER_ICD10Dm HAYFEVER_ICD10Ds HAYFEVER_ICD9 HAYFEVER_SR HAYFEVER_D HAYFEVER_D2 HAYFEVER_GP)
replace HAYFEVER = . if HAYFEVER==1
replace HAYFEVER = 1 if HAYFEVER==2
**********************************************************


**********************************************************
************** HAYFEVER / ASTHMA / ECZEMA ****************
**********************************************************
gen ASTECZRHI=.
replace ASTECZRHI=1 if (ASTHMA==1 | ECZEMA==1 | HAYFEVER==1)
replace ASTECZRHI=0 if (ASTHMA==0 & ECZEMA==0 & HAYFEVER==0)

forval i = 0/3 {
    forval j = 0/2 {
        replace ASTECZRHI = 1 if n_6152_`j'_`i' == 9
    }
}
replace ASTECZRHI = 1 if n_6152_3_0 == 9 | n_6152_3_1 == 9 | n_6152_3_2 == 9

**********************************************************


**********************************************************
************** INFLAMMATORY BOWEL DISEASE ****************
**********************************************************

* ICD10
gen IBD_ICD10 = 0 
forval i = 0/242 {
    forval j = 500/519 {
        replace IBD_ICD10 = 2 if s_41270_0_`i' == "K`j'"
    }
    foreach j in "K523" "K528" {
        replace IBD_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
forval i = 0/242 {
    foreach j in "K529" {
        replace IBD_ICD10 = 1 if IBD_ICD10 == 0 & s_41270_0_`i' == "`j'"
    }
}

**************************

* ICD9
gen IBD_ICD9 = 0 
forval i = 0/46 {
    foreach j in "5550" "5551" "5552" "5559" "5569" {
        replace IBD_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}
forval i = 0/46 {
    foreach j in "558" "5589" "55899" {
        replace IBD_ICD9 = 1 if IBD_ICD9 == 0 & s_41271_0_`i' == "`j'"
    }
}

**************************

* Death register ICD10 - main
gen IBD_ICD10Dm = 0 
forval i = 0/1 {
    forval j = 500/519 {
        replace IBD_ICD10Dm = 2 if s_40001_`i'_0 == "K`j'" 
    }
    foreach j in "K523" "K528" {
        replace IBD_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }    
}
forval i = 0/1 {
    foreach j in "K529" {
        replace IBD_ICD10Dm = 1 if IBD_ICD10Dm == 0 & s_40001_`i'_0 == "`j'"
    }    
}

**************************

* Death register ICD10 - secondary
gen IBD_ICD10Ds = 0 
forval i = 1/14 {
    forval j = 500/519 {
        replace IBD_ICD10Ds = 2 if s_40002_0_`i' == "K`j'" 
    }
    foreach j in "K523" "K528" {
        replace IBD_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }    
}
forval i = 1/9 {
    forval j = 500/519 {
        replace IBD_ICD10Ds = 2 if s_40002_1_`i' == "K`j'" 
    }
    foreach j in "K523" "K528" {
        replace IBD_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }    
}

forval i = 1/14 {
    foreach j in "K529" {
        replace IBD_ICD10Ds = 1 if IBD_ICD10Ds == 0 & s_40002_0_`i' == "`j'"
    }    
}
forval i = 1/9 {
    foreach j in "K529" {
        replace IBD_ICD10Ds = 1 if IBD_ICD10Ds == 0 & s_40002_1_`i' == "`j'"
    }    
}

**************************
egen IBD=rmax(IBD_ICD10 IBD_ICD10Dm IBD_ICD10Ds IBD_ICD9 IBD_GP)
replace IBD = . if IBD==1
replace IBD = 1 if IBD==2
**********************************************************
**********************************************************


**********************************************************
************************* MIGRAINE ***********************
**********************************************************
* In-patient - ICD10
gen MIGRAINE_ICD10=0
forval i = 0/242 {
    foreach j in "G430" "G431" "G432" "G433" "G438" "G439" {
        replace MIGRAINE_ICD10 = 2 if s_41270_0_`i' == "`j'"
    }
}
* Exclude these from controls: headache under "symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"
forval i = 0/242 {
    replace MIGRAINE_ICD10 = 1 if MIGRAINE_ICD10 == 0 & s_41270_0_`i' == "R51"
}

**************************

* Death register ICD10 - main
gen MIGRAINE_ICD10Dm=0
forval i = 0/1 {
    foreach j in "G430" "G431" "G432" "G433" "G438" "G439" {
        replace MIGRAINE_ICD10Dm = 2 if s_40001_`i'_0 == "`j'"
    }
}
* Exclude these from controls: headache under "symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"
forval i = 0/1 {
    replace MIGRAINE_ICD10Dm = 1 if MIGRAINE_ICD10Dm == 0 & s_40001_`i'_0 == "R51"
}

**************************

* Death register ICD10 - secondary
gen MIGRAINE_ICD10Ds=0
forval i = 1/14 {
    foreach j in "G430" "G431" "G432" "G433" "G438" "G439" {
        replace MIGRAINE_ICD10Ds = 2 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "G430" "G431" "G432" "G433" "G438" "G439" {
        replace MIGRAINE_ICD10Ds = 2 if s_40002_1_`i' == "`j'"
    }
}

forval i = 1/14 {
    replace MIGRAINE_ICD10Ds = 1 if MIGRAINE_ICD10Ds == 0 & s_40002_0_`i' == "R51"
}
forval i = 1/9 {
    replace MIGRAINE_ICD10Ds = 1 if MIGRAINE_ICD10Ds == 0 & s_40002_1_`i' == "R51"
}

**************************
* In-patient - ICD9 
gen MIGRAINE_ICD9=0
forval i = 0/46 {
    foreach j in 3460 3461 3462 3468 3469{
        replace MIGRAINE_ICD9 = 2 if s_41271_0_`i' == "`j'"
    }
}
* Exclude these from controls: headache under "symptoms, signs and ill-defined conditions"
forval i = 0/46 {
    replace MIGRAINE_ICD9 = 1 if MIGRAINE_ICD9 == 0 & s_41271_0_`i' == "7840"
}

**************************
gen MIGRAINE_SR=0

* Non-cancer illness self-report
forval i = 0/28 {
    replace MIGRAINE_SR = 2 if (n_20002_0_`i' == 1265)
}
forval i = 0/15 {
    replace MIGRAINE_SR = 2 if (n_20002_1_`i' == 1265 )
}
forval i = 0/33 {
    replace MIGRAINE_SR = 2 if (n_20002_2_`i' == 1265 )
}
forval i = 0/18 {
    replace MIGRAINE_SR = 2 if (n_20002_3_`i' == 1265 )
}
forval i = 0/28 {
    replace MIGRAINE_SR = 1 if MIGRAINE_SR == 0 & n_20002_0_`i' == 1436
}
forval i = 0/15 {
    replace MIGRAINE_SR = 1 if MIGRAINE_SR == 0 & n_20002_1_`i' == 1436
}
forval i = 0/33 {
    replace MIGRAINE_SR = 1 if MIGRAINE_SR == 0 & n_20002_2_`i' == 1436
}
forval i = 0/18 {
    replace MIGRAINE_SR = 1 if MIGRAINE_SR == 0 & n_20002_3_`i' == 1436
}
**************************

* Experience of pain: Have you ever been told by a doctor that you have had migraine?
gen MIGRAINE_P = n_120016_0_0
replace MIGRAINE_P=. if MIGRAINE_P<0
replace MIGRAINE_P=2 if MIGRAINE_P==1   

**************************

egen MIGRAINE=rmax(MIGRAINE_P MIGRAINE_ICD10 MIGRAINE_ICD10Dm MIGRAINE_ICD10Ds MIGRAINE_ICD9 MIGRAINE_SR MIGRAINE_GP)
replace MIGRAINE = . if MIGRAINE==1
replace MIGRAINE = 1 if MIGRAINE==2


**********************************************************



**********************************************************
******************** NEARSIGHTEDNESS  ********************
**********************************************************
/* In-patient - ICD10
gen NEARSIGHTED_ICD10=0
forval i = 0/242 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10 = 1 if s_41270_0_`i' == "`j'"
    }
}

**************************

* Death register ICD10 - main
gen NEARSIGHTED_ICD10Dm=0
forval i = 0/1 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Dm = 1 if s_40001_`i'_0 == "`j'"
    }
}

**************************

* Death register ICD10 - secondary
gen NEARSIGHTED_ICD10Ds=0
forval i = 1/14 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Ds = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Ds = 1 if s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9 
gen NEARSIGHTED_ICD9=0
forval i = 0/46 {
    foreach j in 3671 {
        replace NEARSIGHTED_ICD9 = 1 if s_41271_0_`i' == "`j'"
    }
}

**************************

* Self-report: "Wears glasses or contact lenses" + Why were you prescribed glasses/contacts? 
gen NEARSIGHTED_SR=0
forval i=0/3{
    replace NEARSIGHTED_SR = 1 if n_6147_`i'_0==1
}

**************************

* Daignosis
gen NEARSIGHTED_D=0
replace NEARSIGHTED_D = 1 if n_20262_0_0==1 | n_20262_0_0==2

**************************

replace NEARSIGHTED_GP=1 if NEARSIGHTED_GP==2
egen NEARSIGHTED=rmax(NEARSIGHTED_ICD10 NEARSIGHTED_ICD10Dm NEARSIGHTED_ICD10Ds NEARSIGHTED_ICD9 NEARSIGHTED_SR NEARSIGHTED_D NEARSIGHTED_GP)

*/


**************************
**************************
* New definition
* In-patient - ICD10
gen NEARSIGHTED_ICD10=.
forval i = 0/242 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10 = 1 if s_41270_0_`i' == "`j'"
    }
}

**************************

* Death register ICD10 - main
gen NEARSIGHTED_ICD10Dm=.
forval i = 0/1 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Dm = 1 if s_40001_`i'_0 == "`j'"
    }
}

**************************

* Death register ICD10 - secondary
gen NEARSIGHTED_ICD10Ds=.
forval i = 1/14 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Ds = 1 if s_40002_0_`i' == "`j'"
    }
}
forval i = 1/9 {
    foreach j in "H442" "H521" {
        replace NEARSIGHTED_ICD10Ds = 1 if s_40002_1_`i' == "`j'"
    }
}

**************************

* In-patient - ICD9 
gen NEARSIGHTED_ICD9=.
forval i = 0/46 {
    foreach j in 3671 {
        replace NEARSIGHTED_ICD9 = 1 if s_41271_0_`i' == "`j'"
    }
}

**************************

* Self-report: "Wears glasses or contact lenses" + Why were you prescribed glasses/contacts? 
gen NEARSIGHTED_SR=.
forval i=0/3{
    replace NEARSIGHTED_SR = 1 if n_6147_`i'_0==1
}

**************************

* Daignosis
gen NEARSIGHTED_D=.
replace NEARSIGHTED_D = 1 if n_20262_0_0==1 | n_20262_0_0==2

**************************

replace NEARSIGHTED_GP=1 if NEARSIGHTED_GP==2
replace NEARSIGHTED_GP=. if NEARSIGHTED_GP==0

egen NEARSIGHTED=rmax(NEARSIGHTED_ICD10 NEARSIGHTED_ICD10Dm NEARSIGHTED_ICD10Ds NEARSIGHTED_ICD9 NEARSIGHTED_SR NEARSIGHTED_D NEARSIGHTED_GP) 

forval i=0/3{
    replace NEARSIGHTED=0 if NEARSIGHTED==. & (n_2207_`i'_0==0 | (n_6147_`i'_0!=. & n_6147_`i'_0>1))
}
replace NEARSIGHTED=0 if NEARSIGHTED==. & n_20262_0_0==0

**********************************************************
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
forval i = 0/1 {
    replace PRCA = 1 if (s_40001_`i'_0 == "C61" | s_40001_`i'_0 == "D075") 
}


* Death register ICD10 - secondary
forval i = 1/14 {
    replace PRCA = 1 if (s_40002_0_`i' == "C61" | s_40002_0_`i' == "D075") 
}
forval i = 1/9 {
    replace PRCA = 1 if (s_40002_1_`i' == "C61" | s_40002_1_`i' == "D075") 
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
forval i = 0/21 {
    replace PRCA = 1 if (s_40006_`i'_0 == "C61" | s_40006_`i'_0 == "D075")
}

**************************
* Cancer register ICD9
forval i = 0/8 {
    replace PRCA = 1 if (s_40013_`i'_0 == "1859" | s_40013_`i'_0 == "2334")
}
forval i = 10/12 {
    replace PRCA = 1 if (s_40013_`i'_0 == "1859" | s_40013_`i'_0 == "2334")
}
foreach i in 14 {
    replace PRCA = 1 if (s_40013_`i'_0 == "1859" | s_40013_`i'_0 == "2334")
}

**************************

replace PRCA = . if n_22001_0_0==0
**********************************************************


**********************************************************
********************* TYPE II DIABETES *******************
**********************************************************

* ICD10
gen T2D_ICD10 = 0 
forval i = 0/242 {
    forval j = 110/119 {
        replace T2D_ICD10 = 2 if s_41270_0_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D_ICD10 = 1 if T2D_ICD10 == 0 & s_41270_0_`i' == "E`j'" 
    }
}

* ICD9
gen T2D_ICD9 = 0
forval i = 0/46 {
    foreach j in "25000" "25010" "25020" "25090" {
        replace T2D_ICD9 = 2 if s_41271_0_`i' == "`j'" 
    }
    foreach j in "25009" "25019" "25029" "2503" "2504" "2505" "2506" "2507" "25099" {
        replace T2D_ICD9 = 1 if T2D_ICD9 == 0 & s_41271_0_`i' == "`j'"
    }
}

* Death register ICD10 - main
gen T2D_ICD10Dm = 0
forval i = 0/1 {
    forval j = 110/119 {
        replace T2D_ICD10Dm = 2 if s_40001_`i'_0 == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D_ICD10Dm = 1 if T2D_ICD10Dm == 0 & s_40001_`i'_0 == "E`j'"
    }
}


* Death register ICD10 - secondary
gen T2D_ICD10Ds = 0
forval i = 1/14 {
    forval j = 110/119 {
        replace T2D_ICD10Ds = 2 if s_40002_0_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D_ICD10Ds = 1 if T2D_ICD10Ds == 0 & s_40002_0_`i' == "E`j'"
    }
}
forval i = 1/9 {
    forval j = 110/119 {
        replace T2D_ICD10Ds = 2 if s_40002_1_`i' == "E`j'" 
    }
    forval j = 140/149 {
        replace T2D_ICD10Ds = 1 if T2D_ICD10Ds == 0 & s_40002_1_`i' == "E`j'"
    }
}
**************************
egen T2D=rmax(T2D_ICD10 T2D_ICD10Dm T2D_ICD10Ds T2D_ICD9 T2D_GP)
replace T2D = . if T2D==1
replace T2D = 1 if T2D==2
**********************************************************

**********************************************************

*** SAVE FULL DATASET ***
keep n_eid FID IID Sex* Batch BYEAR partition ALZ ALZnoproxy ALLERGYCAT ALLERGYPOLLEN ALLERGYDUST ASTHMA ASTECZRHI BL_LDL BL_HDL BL_CHOL BL_nonHDL BL_TRYG BPdia BPsys BPpulse BRCA COPD ECZEMA HARDCAD HAYFEVER IBD MIGRAINE NEARSIGHTED PRCA SOFTCAD T2D PC*

save "tmp/pgi_repo_health.dta", replace

**********************************************************************************
**********************************************************************************

****************************************
********* RESIDUALIZE & EXPORT *********
****************************************

local health_pooled_sex NALZ ALZnoproxy ALLERGYCAT ALLERGYPOLLEN ALLERGYDUST ASTHMA ASTECZRHI BL_LDL BL_HDL BL_CHOL BL_nonHDL BL_TRYG BPdia BPsys BPpulse COPD ECZEMA HARDCAD HAYFEVER IBD MIGRAINE NEARSIGHTED SOFTCAD T2D
local health_sex_specific BRCA PRCA
local health ALZ ALZnoproxy ALLERGYCAT ALLERGYPOLLEN ALLERGYDUST ASTHMA ASTECZRHI BL_LDL BL_HDL BL_CHOL BL_nonHDL BL_TRYG BPdia BPsys BPpulse BRCA COPD ECZEMA HARDCAD HAYFEVER IBD MIGRAINE NEARSIGHTED PRCA SOFTCAD T2D

by partition, sort: summarize `health'
/*
foreach partition in 1 2 3 {
    foreach var of varlist `health_pooled_sex' {
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
        }



    ** sex specific phenotypes:
    foreach var of varlist `health_sex_specific' {
        
        qui xi:reg `var' BYEAR* i.Batch PC1-PC40 if partition==`partition'
        predict resid, rstandard
        replace resid=. if partition!=`partition'
        export delimited FID IID resid using "input/UKB_`var'_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }
    
}

foreach partition in 1 2 3 {
    foreach var of varlist ALLERGYCAT ALLERGYPOLLEN ALLERGYDUST ECZEMA{
        
        qui xi:reg `var' BYEAR* Sex* i.Batch PC1-PC40 if partition!=`partition'
        predict resid, rstandard
        replace resid=. if partition==`partition'
        export delimited FID IID resid using "input/UKB_`var'_excl_part`partition'.pheno", noq delim(" ") replace

        drop resid
    }

    

}
*/
*** END ***
log close
