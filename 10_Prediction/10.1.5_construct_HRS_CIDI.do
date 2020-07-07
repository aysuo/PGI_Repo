*----------------------------------------------------------------------------------*
* Constructs CIDI phenotype from HRS data
* Date: 04/09/2018
* Author: Joel Becker, based off original script by Aysu Okbay

* Notes:
*		* for some reason sex variables never found. search around
*----------------------------------------------------------------------------------*


********************************************************
************************ Set-up ************************
********************************************************

clear all
set more off
set maxvar 20000
local inputDataDir = "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS"
local WD = "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction"


********************************************************
******************* Import and merge *******************
********************************************************

import delimited "`inputDataDir'/1995_merged_respondent_level.txt"
save "`WD'/tmp/1995_merged_respondent_level.dta", replace
use "`inputDataDir'/TRK2016TR_R.dta", clear
rename (HHID PN) (hhid pn)
destring hhid pn, replace
merge 1:1 hhid pn using "`WD'/tmp/1995_merged_respondent_level.dta", nogenerate
rename (GENDER d636a d638a d391 d393) (sex_1995 mob_1995 yob_1995 moi_1995 yoi_1995)
keep hhid pn d1006-d1030 d1009-d1017 d1031-d1038 sex_1995 mob_1995 yob_1995 moi_1995 yoi_1995
save "`WD'/tmp/CIDI_1995.dta", replace

clear
import delimited "`inputDataDir'/1996_merged_respondent_level.txt"
save "`WD'/tmp/1996_merged_respondent_level.dta", replace
use "`inputDataDir'/TRK2016TR_R.dta", clear
rename (HHID PN) (hhid pn)
destring hhid pn, replace
merge 1:1 hhid pn using "`WD'/tmp/1996_merged_respondent_level.dta", nogenerate
rename (GENDER e636 e638 e391 e393) (sex_1996 mob_1996 yob_1996 moi_1996 yoi_1996)
keep hhid pn e1006-e1030 e1009-e1018 e1031-e1038 sex_1996 mob_1996 yob_1996 moi_1996 yoi_1996
save "`WD'/tmp/CIDI_1996.dta", replace

clear
import delimited "`inputDataDir'/HRS_1998_data.txt"
rename (f686 f968a f970a f704 f703) (sex_1998 mob_1998 yob_1998 moi_1998 yoi_1998)
keep hhid pn f1323-f1362 sex_1998 mob_1998 yob_1998 moi_1998 yoi_1998
save "`WD'/tmp/CIDI_1998.dta", replace

clear
import delimited "`inputDataDir'/HRS_2000_data.txt"
rename (g757 g1051a g1053a g775 g774) (sex_2000 mob_2000 yob_2000 moi_2000 yoi_2000)
keep hhid pn g1456-g1495 sex_2000 mob_2000 yob_2000 moi_2000 yoi_2000
save "`WD'/tmp/CIDI_2000.dta", replace

clear
import delimited "`inputDataDir'/HRS_2002_data.txt"
rename (hx060_r hx004_r hx067_r ha500 ha501) (sex_2002 mob_2002 yob_2002 moi_2002 yoi_2002)
keep hhid pn hc150-hc182 sex_2002 mob_2002 yob_2002 moi_2002 yoi_2002
save "`WD'/tmp/CIDI_2002.dta", replace

clear
import delimited "`inputDataDir'/HRS_2004_data.txt"
rename (jx060_r jx004_r jx067_r ja500 ja501) (sex_2004 mob_2004 yob_2004 moi_2004 yoi_2004)
keep hhid pn jc150-jc182 sex_2004 mob_2004 yob_2004 moi_2004 yoi_2004
save "`WD'/tmp/CIDI_2004.dta", replace

clear
import delimited "`inputDataDir'/HRS_2006_data.txt"
rename (kx060_r kx004_r kx067_r ka500 ka501) (sex_2006 mob_2006 yob_2006 moi_2006 yoi_2006)
keep hhid pn kc150-kc182 sex_2006 mob_2006 yob_2006 moi_2006 yoi_2006
save "`WD'/tmp/CIDI_2006.dta", replace

clear
import delimited "`inputDataDir'/HRS_2008_data.txt"
rename (lx060_r lx004_r lx067_r la500 la501) (sex_2008 mob_2008 yob_2008 moi_2008 yoi_2008)
keep hhid pn lc150-lc182 sex_2008 mob_2008 yob_2008 moi_2008 yoi_2008
save "`WD'/tmp/CIDI_2008.dta", replace

clear
import delimited "`inputDataDir'/HRS_2010_data.txt"
rename (mx060_r mx004_r mx067_r ma500 ma501) (sex_2010 mob_2010 yob_2010 moi_2010 yoi_2010)
keep hhid pn mc150-mc182 sex_2010 mob_2010 yob_2010 moi_2010 yoi_2010
save "`WD'/tmp/CIDI_2010.dta", replace

clear
import delimited "`inputDataDir'/HRS_2012_data.txt"
rename (nx060_r nx004_r nx067_r na500 na501) (sex_2012 mob_2012 yob_2012 moi_2012 yoi_2012)
keep hhid pn nc150-nc182 sex_2012 mob_2012 yob_2012 moi_2012 yoi_2012
save "`WD'/tmp/CIDI_2012.dta", replace

clear
import delimited "`inputDataDir'/HRS_2014_data.txt"
rename (ox060_r ox004_r ox067_r oa500 oa501) (sex_2014 mob_2014 yob_2014 moi_2014 yoi_2014)
keep hhid pn oc150-oc182 sex_2014 mob_2014 yob_2014 moi_2014 yoi_2014
save "`WD'/tmp/CIDI_2014.dta", replace

clear
import delimited "`inputDataDir'/HRS_2016_data.txt"
rename (px060_r px004_r px067_r pa500 pa501) (sex_2016 mob_2016 yob_2016 moi_2016 yoi_2016)
keep hhid pn pc150-pc182 sex_2016 mob_2016 yob_2016 moi_2016 yoi_2016
save "`WD'/tmp/CIDI_2016.dta", replace

foreach i in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014{
	merge 1:1 hhid pn using "`WD'/tmp/CIDI_`i'.dta", nogenerate
}


********************************************************
******************* Rename variables *******************
********************************************************

rename (d1006 d1028 d1007 d1008 d1029 d1030) (DEPR_screen_1995 ANH_screen_1995 DEPR_intensity_1995 DEPR_freq_1995 ANH_intensity_1995 ANH_freq_1995)
rename (e1006 e1028 e1007 e1008 e1029 e1030) (DEPR_screen_1996 ANH_screen_1996 DEPR_intensity_1996 DEPR_freq_1996 ANH_intensity_1996 ANH_freq_1996)
rename (f1323 f1345 f1324 f1325 f1346 f1347) (DEPR_screen_1998 ANH_screen_1998 DEPR_intensity_1998 DEPR_freq_1998 ANH_intensity_1998 ANH_freq_1998)
rename (g1456 g1478 g1457 g1458 g1479 g1480) (DEPR_screen_2000 ANH_screen_2000 DEPR_intensity_2000 DEPR_freq_2000 ANH_intensity_2000 ANH_freq_2000)

rename (d1009 d1010 d1011 d1012 d1013 d1014 d1015 d1016 d1017) (DEPR_anh_1995 DEPR_tired_1995 DEPR_lossApp_1995 DEPR_incrApp_1995 DEPR_sleep_1995 DEPR_sleepFreq_1995 DEPR_ctrate_1995 DEPR_down_1995 DEPR_death_1995)
rename (e1009 e1010 e1011 e1012 e1013 e1014 e1015 e1016 e1018) (DEPR_anh_1996 DEPR_tired_1996 DEPR_lossApp_1996 DEPR_incrApp_1996 DEPR_sleep_1996 DEPR_sleepFreq_1996 DEPR_ctrate_1996 DEPR_down_1996 DEPR_death_1996)
rename (f1326 f1327 f1328 f1329 f1330 f1331 f1332 f1333 f1334) (DEPR_anh_1998 DEPR_tired_1998 DEPR_lossApp_1998 DEPR_incrApp_1998 DEPR_sleep_1998 DEPR_sleepFreq_1998 DEPR_ctrate_1998 DEPR_down_1998 DEPR_death_1998)
rename (g1459 g1460 g1461 g1462 g1463 g1464 g1465 g1466 g1467) (DEPR_anh_2000 DEPR_tired_2000 DEPR_lossApp_2000 DEPR_incrApp_2000 DEPR_sleep_2000 DEPR_sleepFreq_2000 DEPR_ctrate_2000 DEPR_down_2000 DEPR_death_2000)

rename (d1031 d1032 d1033 d1034 d1035 d1036 d1037 d1038) (ANH_tired_1995 ANH_lossApp_1995 ANH_incrApp_1995 ANH_sleep_1995 ANH_sleepFreq_1995 ANH_ctrate_1995 ANH_down_1995 ANH_death_1995)
rename (e1031 e1032 e1033 e1034 e1035 e1036 e1037 e1038) (ANH_tired_1996 ANH_lossApp_1996 ANH_incrApp_1996 ANH_sleep_1996 ANH_sleepFreq_1996 ANH_ctrate_1996 ANH_down_1996 ANH_death_1996)
rename (f1348 f1349 f1350 f1351 f1352 f1353 f1354 f1355) (ANH_tired_1998 ANH_lossApp_1998 ANH_incrApp_1998 ANH_sleep_1998 ANH_sleepFreq_1998 ANH_ctrate_1998 ANH_down_1998 ANH_death_1998)
rename (g1481 g1482 g1483 g1484 g1485 g1486 g1487 g1488) (ANH_tired_2000 ANH_lossApp_2000 ANH_incrApp_2000 ANH_sleep_2000 ANH_sleepFreq_2000 ANH_ctrate_2000 ANH_down_2000 ANH_death_2000)

local year=2002
foreach i in h j k l m n o p{
	rename `i'c150 DEPR_screen_`year'
	rename `i'c151 DEPR_intensity_`year'
	rename `i'c152 DEPR_freq_`year'

	rename `i'c167 ANH_screen_`year'
	rename `i'c168 ANH_intensity_`year'
	rename `i'c169 ANH_freq_`year'

	rename `i'c153 DEPR_anh_`year'
	rename `i'c154 DEPR_tired_`year'
	rename `i'c155 DEPR_lossApp_`year'
	rename `i'c156 DEPR_incrApp_`year'
	rename `i'c157 DEPR_sleep_`year'
	rename `i'c158 DEPR_sleepFreq_`year'
	rename `i'c159 DEPR_ctrate_`year'
	rename `i'c160 DEPR_down_`year'
	rename `i'c161 DEPR_death_`year'

	rename `i'c170 ANH_tired_`year'
	rename `i'c171 ANH_lossApp_`year'
	rename `i'c172 ANH_incrApp_`year'
	rename `i'c173 ANH_sleep_`year'
	rename `i'c174 ANH_sleepFreq_`year'
	rename `i'c175 ANH_ctrate_`year'
	rename `i'c176 ANH_down_`year'
	rename `i'c177 ANH_death_`year'
	local year=`year'+2
}

gen long hhidpn = 1000*hhid + pn

saveold "`WD'/tmp/CIDI_all.dta", replace version(12)
use "`WD'/tmp/CIDI_all.dta", clear


********************************************************
****************** Construct phenotype *****************
********************************************************

* replace missings
foreach year in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016{
	foreach screen in DEPR ANH{
		foreach var in screen intensity freq tired lossApp incrApp sleep sleepFreq ctrate down death{
			destring `screen'_`var'_`year', replace force
			replace `screen'_`var'_`year' = . if `screen'_`var'_`year' >= 7
		}
	}
	destring DEPR_anh_`year', replace force
	replace DEPR_anh_`year' = . if DEPR_anh_`year' >= 7
}


foreach year in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016{
	* Initiate number of symptoms for individuals that pass either of the screen questions. Anhedonia symptoms start =1 because the screen question itself is considered a symptom question for anhedonia. See http://www.hcp.med.harvard.edu/ncs/ftpdir/cidisf_readme.pdf
	gen symptoms_DEPR_`year'=0
	gen symptoms_ANH_`year'=1
	foreach screen in DEPR ANH{
		* Increment symptom count by 1 for each of the 6 symptoms common among dysphoria and anhedonia.
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_tired_`year'==1
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_lossApp_`year'==1 | `screen'_incrApp_`year'==1
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_sleep_`year'==1 & `screen'_sleepFreq_`year' <=2
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_ctrate_`year'==1
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_down_`year'==1
		replace symptoms_`screen'_`year'=symptoms_`screen'_`year'+1 if `screen'_death_`year'==1
	}
	* Increment dysphoria symptoms by 1 if R has the anhedonia symptom.
	replace symptoms_DEPR_`year'=symptoms_DEPR_`year'+1 if DEPR_anh_`year'==1
}


*Generate CIDI score for each year.
foreach year in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016{
	* Set CIDI to -1 if R answered "No" to both screen questions.
	gen CIDI_`year'=-1  if DEPR_screen_`year'==5 & ANH_screen_`year'==5
	* Set to -1 if R answered "No" to dysphoria, "Yes" to anhedonia, but has failed intensity or frequency requirement of anhedonia.
	replace CIDI_`year'=-1 if DEPR_screen_`year'==5 & ANH_screen_`year'==1 & (ANH_intensity_`year'>2 | ANH_freq_`year'>2)
	* Set to -1 if R answered "Yes" to dysphoria, failed frequency or intensity for dysphoria, and then answered "No" to anhedonia.
	replace CIDI_`year'=-1 if DEPR_screen_`year'==1 & (DEPR_intensity_`year'>2 | DEPR_freq_`year'>2) & ANH_screen_`year'==5
	* Set to -1 if R answered "Yes" to dysphoria, failed frequency or intensity for dysphoria, answered "Yes" to anhedonia, and failed frequency or intensity requirement of anhedonia.
	replace CIDI_`year'=-1 if DEPR_screen_`year'==1 & (DEPR_intensity_`year'>2 | DEPR_freq_`year'>2) & ANH_screen_`year'==1 & (ANH_intensity_`year'>2 | ANH_freq_`year'>2)

	* Set CIDI to number of symptoms for respondents that endorsed one of the screen questions with the appropriate frequency/duration.
	replace CIDI_`year'=symptoms_DEPR_`year' if DEPR_screen_`year'==1 & DEPR_intensity_`year'<=2 & DEPR_freq_`year'<=2
	replace CIDI_`year'=symptoms_ANH_`year' if DEPR_screen_`year'==5 & ANH_screen_`year'==1 & ANH_intensity_`year'<=2 & ANH_freq_`year'<=2
	replace CIDI_`year'=symptoms_ANH_`year' if DEPR_screen_`year'==1 & (DEPR_intensity_`year'>2 | DEPR_freq_`year'>2) & ANH_screen_`year'==1 & ANH_intensity_`year'<=2 & ANH_freq_`year'<=2

	replace CIDI_`year'=. if (DEPR_screen_`year'==1 & DEPR_intensity_`year'<=2 & DEPR_freq_`year'<=2) & (DEPR_anh_`year'==. | DEPR_tired_`year'==. | (DEPR_lossApp_`year'==. & DEPR_incrApp_`year'==.) | DEPR_sleep_`year'==. | (DEPR_sleep_`year'==1 & DEPR_sleepFreq_`year'==.) | DEPR_ctrate_`year'==. | DEPR_down_`year'==. | DEPR_death_`year'==.)
	replace CIDI_`year'=. if (ANH_screen_`year'==1 & ANH_intensity_`year'<=2 & ANH_freq_`year'<=2) & (ANH_tired_`year'==. | (ANH_lossApp_`year'==. & ANH_incrApp_`year'==.) | ANH_sleep_`year'==. | (ANH_sleep_`year'==1 & ANH_sleepFreq_`year'==.) | ANH_ctrate_`year'==. | ANH_down_`year'==. | ANH_death_`year'==.)
}


gen CIDIcount=0
foreach year in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016{
	replace CIDIcount = CIDIcount + 1 if CIDI_`year' !=.
}
drop if CIDIcount==0

keep hhid pn hhidpn CIDI* sex_* mob_* yob_* moi_* yoi_*
saveold "`WD'/tmp/CIDI_all.dta", replace version(12)


********************************************************
************ Merge with crossref and predict ***********
********************************************************

use "`inputDataDir'/HRS_GENOTYPEV2_XREF.dta", clear
destring HHID, replace
destring PN, replace
gen long hhidpn = (1000 * HHID) + PN
merge 1:1 hhidpn using "`WD'/tmp/CIDI_all.dta", keep(match) nogenerate
*merge 1:1 IID using "/var/genetics/dbgap/aokbay/HRS/hwe_eur_keep.txt",  keep(match) nogenerate
* don't need to merge above file because non-eur dropped at PC merge stage

save "`WD'/tmp/CIDISF_eur.dta", replace

*****************************************
clear all
set maxvar 20000
use "`inputDataDir'/randhrs1992_2016v2.dta"
keep hhidpn r*agey_m ragender

rename (r2agey_m r3agey_m r4agey_m r5agey_m r6agey_m r7agey_m r8agey_m r9agey_m r10agey_m r11agey_m) (age_1995 age_1996 age_1998 age_2000 age_2002 age_2004 age_2006 age_2008 age_2010 age_2012)

merge 1:1 hhidpn using "`WD'/tmp/CIDISF_eur.dta", keep(match) nogenerate

* old results:
*  CIDIcount |      Freq.     Percent        Cum.
*------------+-----------------------------------
*          1 |        380        4.41        4.41
*          2 |        883       10.24       14.65
*          3 |      1,090       12.64       27.28
*          4 |      6,267       72.67       99.95
*          5 |          4        0.05      100.00

* new results:
*  CIDIcount |      Freq.     Percent        Cum.
*------------+-----------------------------------
*          1 |        687        4.41        4.41
*          2 |      1,427        9.16       13.57
*          3 |      1,377        8.84       22.42
*          4 |      3,031       19.46       41.88
*          5 |      1,981       12.72       54.60
*          6 |      7,061       45.34       99.94
*          7 |          9        0.06      100.00
*------------+-----------------------------------
*      Total |     15,573      100.00



*********************************************** Diagnosis ************************************************
*Short Form MD Score Probability of CIDI Caseness
*-1 0
*0 0.0001
*1 0.0568
*2 0.2352
*3 0.5542
*4 0.8125
*5 0.8895
*6 0.8895
*7 0.9083

foreach year in 1995 1996 1998 2000 2002 2004 2006 2008 2010 2012{
	gen Depr_prob_`year'=0 if CIDI_`year'==-1
	replace Depr_prob_`year'=0.0001 if CIDI_`year'==0
	replace Depr_prob_`year'=0.0568 if CIDI_`year'==1
	replace Depr_prob_`year'=0.2352 if CIDI_`year'==2
	replace Depr_prob_`year'=0.5542 if CIDI_`year'==3
	replace Depr_prob_`year'=0.8125 if CIDI_`year'==4
	replace Depr_prob_`year'=0.8895 if CIDI_`year'==5
	replace Depr_prob_`year'=0.8895 if CIDI_`year'==6
	replace Depr_prob_`year'=0.9083 if CIDI_`year'==7

	replace age_`year'=. if Depr_prob_`year'==.
}

egen Depr_prob=rmean(Depr_prob*)
egen age=rmean(age_*)
gen age2=age*age

regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==1
predict Z_DEP_1,rstandard
regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==2
predict Z_DEP_2,rstandard
regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==3
predict Z_DEP_3,rstandard
regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==4
predict Z_DEP_4,rstandard
regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==5
predict Z_DEP_5,rstandard
regress Depr_prob age age2 ragender##c.age ragender##c.age2 if CIDIcount==6 | CIDIcount==7
predict Z_DEP_67,rstandard

gen Z_DEP=Z_DEP_1 if CIDIcount==1
replace Z_DEP=Z_DEP_2 if CIDIcount==2
replace Z_DEP=Z_DEP_3 if CIDIcount==3
replace Z_DEP=Z_DEP_4 if CIDIcount==4
replace Z_DEP=Z_DEP_5 if CIDIcount==5
replace Z_DEP=Z_DEP_67 if CIDIcount==6 | CIDIcount==7


********************************************************
**************** Merge with PCs and save ***************
********************************************************

rename LOCAL_ID iid
destring iid, replace
keep iid hhid pn hhidpn Z_*
egen phenotype = std(Z_DEP)
save "`WD'/tmp/HRS_CIDISF_FINAL.dta", replace

export delimited hhidpn phenotype using "`WD'/input/HRS/DEP.pheno", noq delim(",") replace


