*----------------------------------------------------------------------------------*
* Saves WLS data for use in R
* Author: Joel Becker

* Notes:
*   imputed age is ~0.4 greater than actual on average
*   * re-mean to assume surveys were conducted earlier in year?
*----------------------------------------------------------------------------------*


********************************************************
************************ Set-up ************************
********************************************************

clear all
set maxvar 20000
cd "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction/tmp"
use "/disk/genetics3/WLS_DBGAP/phenotype_data/wls_plg_13_06.dta"
_strip_labels *


********************************************************
*************** Rename background columns **************
********************************************************

* survey info
rename (idpriv subject_id rtype) (id_old id_new respondent_type)
destring id_new, generate(id)

* demographic variables
rename (z_sexrsp z_brdxdy inafroa) (sex yob african_american)
replace yob = . if yob < 0
replace sex = . if sex < 0
gen male = 2 - sex
* not excluding self-reported african-americans, will be excluded at scores stage

* can't find date of interview or more detailed age information, so construct using wave year
gen age_57 = 57 - yob if yob != .
gen age_75 = 75 - yob if yob != .
gen age_93 = 93 - yob if yob != .
gen age_04 = 104 - yob if yob != .
gen age_11 = 111 - yob if yob != .
replace age_75 = age_75 + 2 if respondent_type == "s"
replace age_93 = age_93 + 1 if respondent_type == "s"
replace age_04 = age_04 + 1 if respondent_type == "s"

* can now find some age info in documentation
gen age_1957 = agersp
gen age_1975 = z_age75
gen age_1993 = z_ra029re
gen age_2004 = z_ga003re
gen age_2011 = z_ha003re

replace age_1957 = . if age_1957 < 0
replace age_1975 = . if age_1975 < 0
replace age_1993 = . if age_1993 < 0
replace age_2004 = . if age_2004 < 0
replace age_2011 = . if age_2011 < 0

* use derived age variable if given ones don't have value
replace age_1957 = age_57 if age_1957 == .
replace age_1975 = age_75 if age_1975 == .
replace age_1993 = age_93 if age_1993 == .
replace age_2004 = age_04 if age_2004 == .
replace age_2011 = age_11 if age_2011 == .


********************************************************
*************** Rename phenotype columns ***************
********************************************************

* activity
rename (z_mx005rer z_ixe01rer z_jz165rer z_jz168rer) (light_exercise_1993 light_exercise_2004 lightalone_exercise_2011 lighttogether_exercise_2011)
rename (z_mx006rer z_ixe02rer z_jz171rer z_jz174rer) (heavy_exercise_1993 heavy_exercise_2004 heavyalone_exercise_2011 heavytogether_exercise_2011)

* note that heavy_exercise_[1993/2004] are already composite measures, whereas 2011 needs to be created
* scale in 2004-2011 same as 1993, yet range much greater

* adhd
rename z_hh125rek ADHD_2011

* AFB
rename z_agrkd1 AFB_1975
rename z_hd01401 child_dob_2011
*gen AFB_2011 = z_hd01401 - (yob * 12) if z_hd01401 > 0
* first var: participant's age at the birth of the first child, assume in months
* second var: dob of first child

* age first menses
rename (nn001rer z_in190rer) (AFM_1993 AFM_2004)
* 1993 variable is siblings-only so 1994, say 1993 for wrangling reasons later

* agree
rename (z_rh009rec z_mh009rei z_ih009rei z_jh009rei) (agree_phone_1993 agree_1993 agree_2004 agree_2011)
rename (z_rh010rec z_mh010re z_ih010re z_jh010re) (agree_phone_nanswered_1993 agree_nanswered_1993 agree_nanswered_2004 agree_nanswered_2011)
* don't use phone surveys and n-answered later for any phenotypes (see documentation)

* allergies (including hayfever)
rename (z_jx410rer z_jx109rer) (whichallergies_2011 haveallergies_2011)

* asthma
rename (z_mx085rer z_ix085rer z_jx085rer) (asthma_1993 asthma_2004 asthma_2011)

* audit
** commentary on whether Qs included in HRS
* analogy to "In your entire life, have you had at least 12 drinks of any type of alcoholic beverage?" not included
* analogy to "Have you ever taken a drink first thing in the morning to steady your nerves or get rid of a hangover?" not included

* Have you ever drunk alcoholic beverages such as beer, wine, liquor, or mixed alcoholic drinks? Yes/No
rename (z_ru025re z_gu025re z_hu025re) (EverHadDrink_1993 EverHadDrink_2004 EverHadDrink_2011)
* During the last month on how many days did you drink any alcoholic beverages such as beer, wine, liquor, or mixed alcoholic drinks?
rename (z_ru026re z_gu026re z_hu026re) (DaysDrinkingPerWeek_1993 DaysDrinkingPerWeek_2004 DaysDrinkingPerWeek_2011)
*What is the average number of drinks you had on the days you consumed any alcoholic beverages such as beer, wine, liquor, or mixed alcoholic drinks in the past month?
rename (z_ru027re z_gu027re z_hu027re) (DrinksPerDrinkingSession_1993 DrinksPerDrinkingSession_2004 DrinksPerDrinkingSession_2011)
*Number of times you had 5 or more drinks on the same occasion during the last month.
rename (z_ru029re z_gu029re z_hu029re) (DaysHadManyDrinks_1993 DaysHadManyDrinks_2004 DaysHadManyDrinks_2011)
*Have you ever felt bad or guilty about drinking?
rename (z_ru030re z_gu030re z_hu030re) (GuiltyAboutDrinking_1993 GuiltyAboutDrinking_2004 GuiltyAboutDrinking_2011)
* below is analogy to "Have you ever felt that you should cut down on drinking?"
rename (z_ru035re z_gu035re z_hu035re) (WantedCutDrinking_1993 WantedCutDrinking_2004 WantedCutDrinking_2011)
rename (z_ru036re z_gu036re z_hu036re) (WantedCutDrinkingSelf_1993 WantedCutDrinkingSelf_2004 WantedCutDrinkingSelf_2011)
* below is analogy to "Have people ever annoyed you by criticizing your drinking?"
rename (z_ru031re z_gu031re z_hu031re) (DrinkingCriticised_1993 DrinkingCriticised_2004 DrinkingCriticised_2011)
* below have no analogy in HRS
rename (z_ru032re z_gu032re z_hu032re) (DrinkingAnnoyedWork_1993 DrinkingAnnoyedWork_2004 DrinkingAnnoyedWork_2011)
rename (z_ru033re z_gu033re z_hu033re) (DrinkingAnnoyedFamily_1993 DrinkingAnnoyedFamily_2004 DrinkingAnnoyedFamily_2011)

* bmi
rename (srbmi z_mx011rec z_ix011rec z_jx011rec) (bmi_1957 bmi_1993 bmi_2004 bmi_2011)

* consc
rename (z_rh007rec z_mh017rei z_ih017rei z_jh017rei) (consc_phone_1993 consc_1993 consc_2004 consc_2011)
rename (z_rh008rec z_mh018re z_ih018re z_jh018re) (consc_phone_nanswered_1993 consc_nanswered_1993 consc_nanswered_2004 consc_nanswered_2011)

* COPD
rename (z_mx089rer z_ix089rer z_jx089rer) (copd_1993 copd_2004 copd_2011)

* CPD
* items for DID and DO -- clarify this with aysu 
* mx015rer: About how many packs did/do you usually smoke per day then/now?
* ixt08rer: On average, how many packs of cigarettes do you smoke a day?
* jxt08rer: On average, how many packs of cigarettes do you smoke a day?
rename (mx015rer z_ixt08rer z_jxt08rer) (CPD_1993 CPD_2004 CPD_2011)
* condition in smoking now
rename (z_mx013rer z_ix013rec z_jx013rec) (smokenow_1993 smokenow_2004 smokenow_2011)
* unclear if should convert 1993 to cigs, always use packs, or fine to residualise as usual

* depression
rename (z_mu001rec z_iu001rec z_ju001rec) (depr_1993 depr_2004 depr_2011)
rename (z_mu002re z_iu002re z_ju002re) (depr_nanswered_1993 depr_nanswered_2004 depr_nanswered_2011)

* DPW
** variables already renamed in AUDIT

* EA
rename (z_gb103red z_hb103red) (EA_2004 EA_2011)

* eversmoke
rename (z_mx012rer z_ix012rer z_jx012rer) (eversmoke_1993 eversmoke_2004 eversmoke_2011)

* extraversion
rename (z_rh001rec z_mh001rei z_ih001rei z_jh001rei) (extra_phone_1993 extra_1993 extra_2004 extra_2011)
rename (z_rh002rec z_mh002re z_ih002re z_jh002re) (extra_phone_nanswered_1993 extra_nanswered_1993 extra_nanswered_2004 extra_nanswered_2011)

* famsat
* what about jn031rer? aysu said no but bigger sample
* only one var because rest in siblings
rename (z_gb040re) (famsat_2004)

* finsat
rename (z_gp226re z_hp226re) (finsat_2004 finsat_2011)

* friendsat
rename (z_iv201rer z_iv041re) (friendsat1_2004 friendsat2_2004)

* height
rename (z_mx010rec z_ix010rec z_jx010rec) (height_1993 height_2004 height_2011)
*rename (z_mx010rec z_ix010fre z_ix010ire z_ix010rec z_jx010fre z_jx010ire z_jx010rec) (height_1993 height_feet_2004 height_inches_2004 height_2004 height_feet_2011 height_inches_2011 height_2011)
* check concordance (https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=capigrad&module=jmail_health)

* intelligence
rename (ghnrs_bm shnrs_t) (intelligence_1957 intelligence_1975)
* note intelligence measure actually from 1977

* left out social
rename (z_sofrno z_mz023rer z_iz023rer z_jz023rer) (withfriends_1975 withfriends_1993 withfriends_2004 withfriends_2011)
rename (z_mn032rer z_in032rer z_jn032rer) (fewfriends_1993 fewfriends_2004 fewfriends_2011)
rename (z_mn033rer z_in033rer z_jn033rer) (otherpeople_1993 otherpeople_2004 otherpeople_2011)

* lonely
rename (z_mu008rer z_iu008rer z_ju008rer) (lonely_1993 lonely_2004 lonely_2011)

* migraine
* have phenotype for "How often have you had headaches in the past 6 months?" (e.g. jx026rer)
rename (mx026rer z_ix026rer z_jx026rer) (migraine_1993 migraine_2004 migraine_2011)
* as well as "Have you had headaches in the past 6 months?" (e.g. jx025rer)
* all here https://www.ssc.wisc.edu/wlsresearch/documentation/browse/index.php?label=headache&variable=&wave_103_105=on&wave_106=on&wave_108=on&wave_112_113=on&searchButton=Search

* nearsighted
* said no because confusing Qs and not enough data
* e.g. check here: https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=capigrad&module=jhealth (hx203lre)

* neb
rename (z_hd301kd z_hd001kd z_hd201kd) (NEB_since75 NEB_11CAPI NEB_CAPIMOSAQ)

* neuroticism
rename (z_rh005rec z_mh025rei z_ih025rei z_jh025rei) (neur_phone_1993 neur_1993 neur_2004 neur_2011)
rename (z_rh006rec z_mh026re z_ih026re z_jh026re) (neur_phone_nanswered_1993 neur_nanswered_1993 neur_nanswered_2004 neur_nanswered_2011)

* openness
rename (z_rh003rec z_mh032rei z_ih032rei z_jh032rei) (open_phone_1993 open_1993 open_2004 open_2011)
rename (z_rh004rec z_mh033re z_ih033re z_jh033re) (open_phone_nanswered_1993 open_nanswered_1993 open_nanswered_2004 open_nanswered_2011)

* religiosity
rename (z_bkxrl3 z_rl004rec z_gl004rec z_hl004rec z_jl004rec) (relig_1975 relig_1993 relig_2004 relig_mail_2011 relig_phone_2011)
* for 1993 not sure if rl004rec or rl002rec + rl003red (think effectively the same)
* for 2011 not sure if hl004rec or jl004rec (come from different surveys)
* gen rel_2011 = (z_hl004rec >= 0 & z_jl004rec >= 0)

* risk

* self-rated health
rename (z_mx001rer z_ix001rer z_jx001rer) (selfhealth_1993 selfhealth_2004 selfhealth_2011)

* SWB
rename (z_mu006rer z_iu006rer z_ju006rer z_mu009rer z_iu009rer z_ju009rer) (happy_1993 happy_2004 happy_2011 enlife_1993 enlife_2004 enlife_2011)

* worksat
rename (z_jbcrsa z_rg044jjc z_gg044jjc z_hg044jjc) (worksat_1975 worksat_1993 worksat_2004 worksat_2011)


********************************************************
************** Recode phenotype variables **************
********************************************************

* activity
replace light_exercise_1993 = . if light_exercise_1993 < 0
recode light_exercise_1993 (1 = 15) (2 = 6) (3 = 2) (4 = 0.5)

replace heavy_exercise_1993 = . if heavy_exercise_1993 < 0
recode heavy_exercise_1993 (1 = 15) (2 = 6) (3 = 2) (4 = 0.5)

replace light_exercise_2004 = . if light_exercise_2004 < 0
replace heavy_exercise_2004 = . if heavy_exercise_2004 < 0

gen light_exercise_2011 = .
replace lightalone_exercise_2011 = . if lightalone_exercise_2011 < 0
replace lighttogether_exercise_2011 = . if lighttogether_exercise_2011 < 0
replace light_exercise_2011 = (lightalone_exercise_2011 + lighttogether_exercise_2011)

gen heavy_exercise_2011 = .
replace heavyalone_exercise_2011 = . if heavyalone_exercise_2011 < 0
replace heavytogether_exercise_2011 = . if heavytogether_exercise_2011 < 0
replace heavy_exercise_2011 = (heavyalone_exercise_2011 + heavytogether_exercise_2011)

gen activity_1993 = .
replace activity_1993 = (2 * light_exercise_1993) + (8 * heavy_exercise_1993)

gen activity_2004 = .
replace activity_2004 = (2 * light_exercise_2004) + (8 * heavy_exercise_2004)

gen activity_2011 = .
replace activity_2011 = (2 * light_exercise_2011) + (8 * heavy_exercise_2011)
* assume throughout that exercising once = 1 hour

* adhd
replace ADHD_2011 = . if ADHD_2011 < 0
replace ADHD_2011 = 2 - ADHD_2011
* sample is tiny! 256

* AFB
replace AFB_1975 = . if AFB_1975 < 0
replace AFB_1975 = AFB_1975 / 10
* Note: Implied decimal between second and third digits (e.g. 15.0).
replace child_dob_2011 = . if child_dob_2011 < 0
replace child_dob_2011 = 1900 + (child_dob_2011 / 12)
/* CENTURY MONTH (YEAR in public version) */
gen AFB_2011 = child_dob_2011 - (1900 + yob + 0.5)
replace AFB_2011 = . if AFB_2011 < 0
* temporary fix

*replace AFB_2011 = 6 + (AFB_2011 / 12)
* AFB_2011 is given as 'CENTURY MONTH' (presumably since jan 1900, would give 1933 to 2004)
* if 1975 same, would be 10 to 30
* 1975 has implicit decimal between 3rd and 4th digit!

* age first menses
replace AFM_1993 = . if AFM_1993 < 0
replace AFM_2004 = . if AFM_2004 < 0
replace AFM_2004 = . if AFM_2004 < 5 | AFM_2004 > 30

* agree
replace agree_1993 = . if agree_1993 < 0
replace agree_2004 = . if agree_2004 < 0
replace agree_2011 = . if agree_2011 < 0

* allergies
gen allergy_2011 = 0
replace allergy_2011 = . if haveallergies_2011 < 0
replace allergy_2011 = . if whichallergies_2011 < 0 | whichallergies_2011 == 22

gen anyallergy_2011 = .
replace anyallergy_2011 = 0 if haveallergies_2011 == 2
replace anyallergy_2011 = 1 if haveallergies_2011 == 1

gen pollen_2011 = allergy_2011
replace pollen_2011 = 1 if whichallergies_2011 == 3

gen cat_2011 = allergy_2011
replace cat_2011 = 1 if whichallergies_2011 == 20 | whichallergies_2011 == 2

gen hayfever_2011 = allergy_2011
replace hayfever_2011 = 1 if whichallergies_2011 == 14

gen dust_2011 = allergy_2011
replace dust_2011 = 1 if whichallergies_2011 == 6


* asthma
replace asthma_1993 = 2 - asthma_1993
replace asthma_1993 = . if asthma_1993 > 1

replace asthma_2004 = 2 - asthma_2004
replace asthma_2004 = . if asthma_2004 > 1

replace asthma_2011 = 2 - asthma_2011
replace asthma_2011 = . if asthma_2011 > 1

* asthmahayfever
gen asthmahayfever_2011 = .
replace asthmahayfever_2011 = 0 if asthma_2011 == 0 & hayfever_2011 == 0
replace asthmahayfever_2011 = 1 if asthma_2011 == 1 | hayfever_2011 == 1

* audit
** note, use same thresholds as for HRS construction
foreach i in 1993 2004 2011 {
  * housekeeping
  display `i'

  * discretise binary variables
  replace EverHadDrink_`i' = . if EverHadDrink_`i' < 0
  replace EverHadDrink_`i' = 2 - EverHadDrink_`i'

  replace GuiltyAboutDrinking_`i' = . if GuiltyAboutDrinking_`i' < 0
  replace GuiltyAboutDrinking_`i' = 2 - GuiltyAboutDrinking_`i'

  replace WantedCutDrinking_`i' = . if WantedCutDrinking_`i' < 0
  replace WantedCutDrinking_`i' = 2 - WantedCutDrinking_`i'
  replace WantedCutDrinkingSelf_`i' = . if WantedCutDrinkingSelf_`i' < 0
  recode WantedCutDrinkingSelf_`i' (1 3 = 1) (2 = 0)
  replace WantedCutDrinking_`i' = 0 if WantedCutDrinkingSelf_`i' == 0

  replace DrinkingCriticised_`i' = . if DrinkingCriticised_`i' < 0
  replace DrinkingCriticised_`i' = 2 - DrinkingCriticised_`i'

  replace DrinkingAnnoyedWork_`i' = . if DrinkingAnnoyedWork_`i' < 0
  replace DrinkingAnnoyedWork_`i' = 2 - DrinkingAnnoyedWork_`i'

  replace DrinkingAnnoyedFamily_`i' = . if DrinkingAnnoyedFamily_`i' < 0
  replace DrinkingAnnoyedFamily_`i' = 2 - DrinkingAnnoyedFamily_`i'

  * discretise quasi-continuous variables using cut-offs from HRS
  ** note, monthly variable
  replace DaysDrinkingPerWeek_`i' = . if DaysDrinkingPerWeek_`i' < 0
  replace DaysDrinkingPerWeek_`i' = (DaysDrinkingPerWeek_`i' / 30) * 7
  gen DaysDrinkingPerWeek_binary_`i' = DaysDrinkingPerWeek_`i'
  replace DaysDrinkingPerWeek_binary_`i' = 0 if DaysDrinkingPerWeek_binary_`i' <= 1
  replace DaysDrinkingPerWeek_binary_`i' = 1 if DaysDrinkingPerWeek_binary_`i' > 1

  ** note, daily variable
  replace DrinksPerDrinkingSession_`i' = . if DrinksPerDrinkingSession_`i' < 0
  gen DrinksPerSession_binary_`i' = DrinksPerDrinkingSession_`i'
  replace DrinksPerSession_binary_`i' = 0 if DrinksPerSession_binary_`i' <= 4
  replace DrinksPerSession_binary_`i' = 1 if DrinksPerSession_binary_`i' > 4

  ** note, monthly variable
  replace DaysHadManyDrinks_`i' = . if DaysHadManyDrinks_`i' < 0
  replace DaysHadManyDrinks_`i' = (DaysHadManyDrinks_`i' * 3)
  replace DaysHadManyDrinks_`i' = 0 if DaysHadManyDrinks_`i' <= 3
  replace DaysHadManyDrinks_`i' = 1 if DaysHadManyDrinks_`i' > 3

  * sum values and empty values across variables
  gen number_na_audit_`i' = 0
  gen audit_`i' = 0

  foreach j in GuiltyAboutDrinking WantedCutDrinkingSelf WantedCutDrinking DrinkingCriticised DrinkingAnnoyedWork DrinkingAnnoyedFamily DaysDrinkingPerWeek_binary DrinksPerSession_binary DaysHadManyDrinks {
    replace number_na_audit_`i' = number_na_audit_`i' + (`j'_`i' == .)
    replace audit_`i' = audit_`i' + `j'_`i' if `j'_`i' != .
  }

  * set to missing if less than minimum number of variables
  replace audit_`i' = . if number_na_audit_`i' > 4

  * set to zero if filter question set to zero
  replace audit_`i' = 0 if EverHadDrink_`i' == 0
}

* bmi
replace bmi_1993 = . if bmi_1993 < 0
replace bmi_2004 = . if bmi_2004 < 0
replace bmi_2011 = . if bmi_2011 < 0
* note 1957 scale is weird

* consc
replace consc_1993 = . if consc_1993 < 0
replace consc_2004 = . if consc_2004 < 0
replace consc_2011 = . if consc_2011 < 0
* want to do something with n_answered? espesh since numbers so high, maybe need to average/check documentation

* COPD
replace copd_1993 = . if copd_1993 < 0
replace copd_1993 = 2 - copd_1993

replace copd_2004 = . if copd_2004 < 0
replace copd_2004 = 2 - copd_2004

replace copd_2011 = . if copd_2011 < 0
replace copd_2011 = 2 - copd_2011

* CPD
replace CPD_1993 = . if CPD_1993 < 0
replace CPD_1993 = . if CPD_1993 > 4
* excluding pipes above, ask aysu!
*replace smokenow_1993 = . if smokenow_1993 < 0
*replace smokenow_1993 = 2 - smokenow_1993
*replace CPD_1993 = . if smokenow_1993 == 0

replace CPD_2004 = . if CPD_2004 < 0
replace CPD_2004 = . if CPD_2004 > 10
*replace smokenow_2004 = . if smokenow_2004 < 0
*replace smokenow_2004 = 2 - smokenow_2004
*replace CPD_2004 = . if smokenow_2004 == 0

replace CPD_2011 = . if CPD_2011 < 0
replace CPD_2011 = . if CPD_2011 > 10
*replace smokenow_2011 = . if smokenow_2011 < 0
*replace smokenow_2011 = 2 - smokenow_2011
*replace CPD_2011 = . if smokenow_2011 == 0

* depression
replace depr_1993 = . if depr_1993 < 0
replace depr_2004 = . if depr_2004 < 0
replace depr_2011 = . if depr_2011 < 0

* dpw
gen dpw_1993 = DaysDrinkingPerWeek_1993 * DrinksPerDrinkingSession_1993
gen dpw_2004 = DaysDrinkingPerWeek_2004 * DrinksPerDrinkingSession_2004
gen dpw_2011 = DaysDrinkingPerWeek_2011 * DrinksPerDrinkingSession_2011

* EA
replace EA_2004 = . if EA_2004 < 0
replace EA_2011 = . if EA_2011 < 0
replace EA_2004 = . if EA_2004 >= 21
replace EA_2011 = . if EA_2011 >= 21

* eversmoke
replace eversmoke_1993 = . if eversmoke_1993 < 0
replace eversmoke_1993 = 2 - eversmoke_1993
replace eversmoke_2004 = . if eversmoke_2004 < 0
replace eversmoke_2004 = 2 - eversmoke_2004
replace eversmoke_2011 = . if eversmoke_2011 < 0
replace eversmoke_2011 = 2 - eversmoke_2011

* extraversion
replace extra_1993 = . if extra_1993 < 0
replace extra_2004 = . if extra_2004 < 0
replace extra_2011 = . if extra_2011 < 0

* famsat
replace famsat_2004 = . if famsat_2004 < 0

* finsat
replace finsat_2004 = . if finsat_2004 < 0
replace finsat_2004 = 5 - finsat_2004
replace finsat_2011 = . if finsat_2011 < 0
replace finsat_2011 = 5 - finsat_2011

* friendsat
replace friendsat1_2004 = . if friendsat1_2004 < 0

replace friendsat2_2004 = . if friendsat2_2004 < 0
replace friendsat2_2004 = 2 - friendsat2_2004

* height
replace height_1993 = . if height_1993 < 0
replace height_1993 = height_1993 / 100
* above now in inches
replace height_2004 = . if height_2004 < 0
replace height_2011 = . if height_2011 < 0

* intelligence
replace intelligence_1957 = . if intelligence_1957 < 0
replace intelligence_1975 = . if intelligence_1975 < 0
* note intelligence measure actually from 1977

* leftoutsocial
gen leftoutsocial_count_1993 = 0
gen leftoutsocial_count_2004 = 0
gen leftoutsocial_count_2011 = 0

replace withfriends_1993 = . if withfriends_1993 < 0
replace withfriends_2004 = . if withfriends_2004 < 0
replace withfriends_2011 = . if withfriends_2011 < 0

replace fewfriends_1993 = . if fewfriends_1993 < 0
replace fewfriends_2004 = . if fewfriends_2004 < 0
replace fewfriends_2011 = . if fewfriends_2011 < 0

replace otherpeople_1993 = . if otherpeople_1993 < 0
replace otherpeople_2004 = . if otherpeople_2004 < 0
replace otherpeople_2011 = . if otherpeople_2011 < 0

quietly: sum withfriends_1993
local mean = r(mean)
local sd = r(sd)
replace withfriends_1993 = (withfriends_1993 - `mean') / `sd'
quietly: sum withfriends_2004
local mean = r(mean)
local sd = r(sd)
replace withfriends_2004 = (withfriends_2004 - `mean') / `sd'
quietly: sum withfriends_2011
local mean = r(mean)
local sd = r(sd)
replace withfriends_2011 = (withfriends_2011 - `mean') / `sd'

quietly: sum fewfriends_1993
local mean = r(mean)
local sd = r(sd)
replace fewfriends_1993 = (fewfriends_1993 - `mean') / `sd'
quietly: sum fewfriends_2004
local mean = r(mean)
local sd = r(sd)
replace fewfriends_2004 = (fewfriends_2004 - `mean') / `sd'
quietly: sum fewfriends_2011
local mean = r(mean)
local sd = r(sd)
replace fewfriends_2011 = (fewfriends_2011 - `mean') / `sd'

quietly: sum otherpeople_1993
local mean = r(mean)
local sd = r(sd)
replace otherpeople_1993 = (otherpeople_1993 - `mean') / `sd'
quietly: sum otherpeople_2004
local mean = r(mean)
local sd = r(sd)
replace otherpeople_2004 = (otherpeople_2004 - `mean') / `sd'
quietly: sum otherpeople_2011
local mean = r(mean)
local sd = r(sd)
replace otherpeople_2011 = (otherpeople_2011 - `mean') / `sd'

replace leftoutsocial_count_1993 = leftoutsocial_count_1993 + 1 if withfriends_1993 != .
replace leftoutsocial_count_2004 = leftoutsocial_count_2004 + 1 if withfriends_2004 != .
replace leftoutsocial_count_2011 = leftoutsocial_count_2011 + 1 if withfriends_2011 != .

replace leftoutsocial_count_1993 = leftoutsocial_count_1993 + 1 if fewfriends_1993 != .
replace leftoutsocial_count_2004 = leftoutsocial_count_2004 + 1 if fewfriends_2004 != .
replace leftoutsocial_count_2011 = leftoutsocial_count_2011 + 1 if fewfriends_2011 != .

replace leftoutsocial_count_1993 = leftoutsocial_count_1993 + 1 if otherpeople_1993 != .
replace leftoutsocial_count_2004 = leftoutsocial_count_2004 + 1 if otherpeople_2004 != .
replace leftoutsocial_count_2011 = leftoutsocial_count_2011 + 1 if otherpeople_2011 != .

egen leftoutsocial_1993 = rmean(withfriends_1993 fewfriends_1993 otherpeople_1993) if leftoutsocial_count_1993>1
egen leftoutsocial_2004 = rmean(withfriends_2004 fewfriends_2004 otherpeople_2004) if leftoutsocial_count_2004>1
egen leftoutsocial_2011 = rmean(withfriends_2011 fewfriends_2011 otherpeople_2011) if leftoutsocial_count_2011>1

*foreach year in 1993 2004 2011 {
*  gen leftoutsocial_count_`year' = 0
*  foreach var in withfriends fewfriends otherpeople {
*    replace `var'_`year' = . if `var'_`year' < 0
*
*    quietly: sum `var'_`year'
*    local mean = r(mean)
*    local sd = r(sd)
*    replace `var'_`year' = (`var'_`year' - `mean') / `sd'
*
*    replace leftoutsocial_count_`year' = leftoutsocial_count_`year' + 1 if `var' != .
*  }
*  egen leftoutsocial_`year' = rmean(withfriends_`year' fewfriends_`year' otherpeople_`year') if leftoutsocial_count>1
*}

* lonely
replace lonely_1993 = . if lonely_1993 < 0
replace lonely_2004 = . if lonely_2004 < 0
replace lonely_2011 = . if lonely_2011 < 0

* migraine
replace migraine_1993 = 0 if mx025rer == 2
replace migraine_1993 = . if migraine_1993 < 0
replace migraine_2004 = migraine_2004 - 1
replace migraine_2004 = . if migraine_2004 < 0
replace migraine_2011 = migraine_2011 - 1
replace migraine_2011 = . if migraine_2011 < 0

* NEB
gen NEB = NEB_since75

* neuroticism
replace neur_1993 = . if neur_1993 < 0
*replace neur_phone_1993 = . if neur_phone_1993 < 0
replace neur_2004 = . if neur_2004 < 0
replace neur_2011 = . if neur_2011 < 0
* note for (e.g.) mh025rei says:
* Note: MH025REI is coded with a sum if at least three of the five component items recieved a valid response. Missing responses
*       are imputed to the mean of remaining valid items prior to summing. See MH026RE for the number of valid responses.
* thus no need to alter missings or normalise -- seems already done
* besides, see (e.g.) mh026re, almost all answered all Qs

* openness
replace open_1993 = . if open_1993 < 0
replace open_2004 = . if open_2004 < 0
replace open_2011 = . if open_2011 < 0

* religiosity
replace relig_1975 = . if relig_1975 < 0
replace relig_1993 = . if relig_1993 < 0
replace relig_2004 = . if relig_2004 < 0
replace relig_mail_2011 = . if relig_mail_2011 < 0
replace relig_phone_2011 = . if relig_phone_2011 < 0
egen relig_2011 = rmean(relig_mail_2011 relig_phone_2011)

* risk
** gain domain
forvalues i = 2(1)8 {
  * housekeeping
  display `i'
  rename z_jstk0`i're realstakes5_`i'_2011
  * remove NAs
  replace realstakes5_`i'_2011 = . if realstakes5_`i'_2011 < 0
  * recode so {0 = certain, 1 = risky}
  replace realstakes5_`i'_2011 = realstakes5_`i'_2011 - 1
}
rename z_jstk09re realstakes9_9_2011
replace realstakes9_9_2011 = . if realstakes9_9_2011 < 0
replace realstakes9_9_2011 = realstakes9_9_2011 - 1
forvalues i = 10(1)15 {
  * housekeeping
  display `i'
  rename z_jstk`i're realstakes9_`i'_2011
  * remove NAs
  replace realstakes9_`i'_2011 = . if realstakes9_`i'_2011 < 0
  * recode so {0 = certain, 1 = risky}
  replace realstakes9_`i'_2011 = realstakes9_`i'_2011 - 1
}
forvalues i = 16(1)22 {
  * housekeeping
  display `i'
  rename z_jstk`i're realstakes11_`i'_2011
  * remove NAs
  replace realstakes11_`i'_2011 = . if realstakes11_`i'_2011 < 0
  * recode so {0 = certain, 1 = risky}
  replace realstakes11_`i'_2011 = realstakes11_`i'_2011 - 1
}

egen risk5_2011 = rmean(realstakes5_*_2011)
egen risk9_2011 = rmean(realstakes9_*_2011)
egen risk11_2011 = rmean(realstakes11_*_2011)

** loss domain
forvalues i = 25(1)31 {
  * housekeeping
  display `i'
  rename z_jstk`i're realstakes_losing5_`i'_2011
  * remove NAs
  replace realstakes_losing5_`i'_2011 = . if realstakes_losing5_`i'_2011 < 0
  * recode so {1 = certain, 0 = risky}
  replace realstakes_losing5_`i'_2011 = realstakes_losing5_`i'_2011 - 1
}
forvalues i = 32(1)38 {
  * housekeeping
  display `i'
  rename z_jstk`i're realstakes_losing9_`i'_2011
  * remove NAs
  replace realstakes_losing9_`i'_2011 = . if realstakes_losing9_`i'_2011 < 0
  * recode so {1 = certain, 0 = risky}
  replace realstakes_losing9_`i'_2011 = realstakes_losing9_`i'_2011 - 1
}
forvalues i = 39(1)45 {
  * housekeeping
  display `i'
  rename z_jstk`i're realstakes_losing11_`i'_2011
  * remove NAs
  replace realstakes_losing11_`i'_2011 = . if realstakes_losing11_`i'_2011 < 0
  * recode so {1 = certain, 0 = risky}
  replace realstakes_losing11_`i'_2011 = realstakes_losing11_`i'_2011 - 1
}

egen risklosing5_2011 = rmean(realstakes_losing5_*_2011)
egen risklosing9_2011 = rmean(realstakes_losing9_*_2011)
egen risklosing11_2011 = rmean(realstakes_losing11_*_2011)

* selfhealth
replace selfhealth_1993 = . if selfhealth_1993 < 0
replace selfhealth_2004 = . if selfhealth_2004 < 0
replace selfhealth_2011 = . if selfhealth_2011 < 0

* SWB
replace happy_1993 = . if happy_1993 < 0
replace happy_2004 = . if happy_2004 < 0
replace happy_2011 = . if happy_2011 < 0
replace enlife_1993 = . if enlife_1993 < 0
replace enlife_2004 = . if enlife_2004 < 0
replace enlife_2011 = . if enlife_2011 < 0
egen SWB_1993 = rmean(happy_1993 enlife_1993)
egen SWB_2004 = rmean(happy_2004 enlife_2004)
egen SWB_2011 = rmean(happy_2011 enlife_2011)

* worksat
replace worksat_1975 = . if worksat_1975 < 0

replace worksat_1993 = . if worksat_1993 < 0
replace worksat_1993 = 5 - worksat_1993

replace worksat_2004 = . if worksat_2004 < 0
replace worksat_2004 = 5 - worksat_2004

replace worksat_2011 = . if worksat_2011 < 0
replace worksat_2011 = 5 - worksat_2011


********************************************************
*********************** Save data **********************
********************************************************

keep id id_old respondent_type male yob african_american NEB *_19* *_20* age*
outsheet using "WLS_renamed.csv", replace comma nolabel
saveold "WLS_renamed.dta", replace version(12)
