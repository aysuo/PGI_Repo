*----------------------------------------------------------------------------------*
* Constructs SWB phenotype from HRS data
* Date: 04/16/2018
* Author: Joel Becker, based off original script by Aysu Okbay

* Notes:
*		* also 2004 data for some LS questions
*----------------------------------------------------------------------------------*


********************************************************
************************ Set-up ************************
********************************************************

clear all
set more off
set maxvar 20000

local inputDataDir="`1'"

********************************************************
******************* Import and merge *******************
********************************************************

clear
import delimited "`inputDataDir'/HRS_2006_data.txt"
rename (kx060_r kx004_r kx067_r ka500 ka501) (sex_2006 mob_2006 yob_2006 moi_2006 yoi_2006)
keep hhid pn klb027* klb003* sex_2006 mob_2006 yob_2006 moi_2006 yoi_2006
save "tmp/SWB_2006.dta", replace

clear
import delimited "`inputDataDir'/HRS_2008_data.txt"
rename (lx060_r lx004_r lx067_r la500 la501) (sex_2008 mob_2008 yob_2008 moi_2008 yoi_2008)
keep hhid pn llb027* llb003* sex_2008 mob_2008 yob_2008 moi_2008 yoi_2008
save "tmp/SWB_2008.dta", replace

clear
import delimited "`inputDataDir'/HRS_2010_data.txt"
rename (mx060_r mx004_r mx067_r ma500 ma501) (sex_2010 mob_2010 yob_2010 moi_2010 yoi_2010)
keep hhid pn mlb027* mlb003* sex_2010 mob_2010 yob_2010 moi_2010 yoi_2010
save "tmp/SWB_2010.dta", replace

clear
import delimited "`inputDataDir'/HRS_2012_data.txt"
rename (nx060_r nx004_r nx067_r na500 na501) (sex_2012 mob_2012 yob_2012 moi_2012 yoi_2012)
keep hhid pn nlb027* nlb003* sex_2012 mob_2012 yob_2012 moi_2012 yoi_2012
save "tmp/SWB_2012.dta", replace

clear
import delimited "`inputDataDir'/HRS_2014_data.txt"
rename (ox060_r ox004_r ox067_r oa500 oa501) (sex_2014 mob_2014 yob_2014 moi_2014 yoi_2014)
keep hhid pn olb026* olb002* sex_2014 mob_2014 yob_2014 moi_2014 yoi_2014
save "tmp/SWB_2014.dta", replace

clear
import delimited "`inputDataDir'/HRS_2016_data.txt"
rename (px060_r px004_r px067_r pa500 pa501) (sex_2016 mob_2016 yob_2016 moi_2016 yoi_2016)
keep hhid pn plb026* plb002* sex_2016 mob_2016 yob_2016 moi_2016 yoi_2016
save "tmp/SWB_2016.dta", replace

foreach i in 2006 2008 2010 2012 2014{
	merge 1:1 hhid pn using "tmp/SWB_`i'.dta", nogenerate
}

gen long hhidpn = 1000*hhid + pn

saveold "tmp/SWB_all.dta", replace version(12)
use "tmp/SWB_all.dta", clear


********************************************************
****************** Construct phenotype *****************
********************************************************

***** 2006
*** PA
gen missingPA=0
foreach i in klb027a klb027b klb027c klb027d klb027e klb027f klb027g klb027h {
	destring `i', replace force
	replace missingPA=missingPA+1 if `i'==.
}

recode klb027a klb027b klb027c klb027d klb027e klb027f klb027g klb027h (1 = 5) (2 = 4) (5 = 1) (4 = 2)
egen KL_PA=rmean(klb027a klb027b klb027c klb027d klb027e klb027f klb027g klb027h) if missingPA < 5

*** LS
gen missingLS=0
foreach i in klb003a klb003b klb003c klb003d klb003e {
	destring `i', replace force
	replace missingLS=missingLS+1 if `i'==.
}

egen KL_LS=rmean(klb003a klb003b klb003c klb003d klb003e) if missingLS < 3
drop missing*


***** 2008-2012
foreach i in ll ml nl{
	*** PA
	gen missingPA=0

	foreach j in `i'b027c `i'b027d `i'b027f `i'b027g `i'b027h `i'b027k `i'b027p `i'b027q `i'b027t `i'b027u `i'b027v `i'b027x `i'b027y {
		destring `j', replace force
		replace missingPA=missingPA+1 if `j'==.
	}

	recode `i'b027c `i'b027d `i'b027f `i'b027g `i'b027h `i'b027k `i'b027p `i'b027q `i'b027t `i'b027u `i'b027v `i'b027x `i'b027y (1 = 5) (2 = 4) (5 = 1) (4 = 2)
	egen `i'_PA=rmean(`i'b027c `i'b027d `i'b027f `i'b027g `i'b027h `i'b027k `i'b027p `i'b027q `i'b027t `i'b027u `i'b027v `i'b027x `i'b027y) if missingPA < 7

	*** LS
	gen missingLS=0

	foreach j in `i'b003a `i'b003b `i'b003c `i'b003d `i'b003e {
		destring `j', replace force
		replace missingLS=missingLS+1 if `j'==.
	}

	egen `i'_LS=rmean(`i'b003a `i'b003b `i'b003c `i'b003d `i'b003e) if missingLS < 3
	drop missing*
}


***** 2014-2016
foreach i in ol pl{
	*** PA
	gen missingPA=0

	foreach j in `i'b026c `i'b026d `i'b026f `i'b026g `i'b026h `i'b026k `i'b026p `i'b026q `i'b026t `i'b026u `i'b026v `i'b026x `i'b026y {
		destring `j', replace force
		replace missingPA=missingPA+1 if `j'==.
	}

	recode `i'b026c `i'b026d `i'b026f `i'b026g `i'b026h `i'b026k `i'b026p `i'b026q `i'b026t `i'b026u `i'b026v `i'b026x `i'b026y (1 = 5) (2 = 4) (5 = 1) (4 = 2)
	egen `i'_PA=rmean(`i'b026c `i'b026d `i'b026f `i'b026g `i'b026h `i'b026k `i'b026p `i'b026q `i'b026t `i'b026u `i'b026v `i'b026x `i'b026y) if missingPA < 7

	*** LS
	gen missingLS=0

	foreach j in `i'b002a `i'b002b `i'b002c `i'b002d `i'b002e {
		destring `j', replace force
		replace missingLS=missingLS+1 if `j'==.
	}

	egen `i'_LS=rmean(`i'b002a `i'b002b `i'b002c `i'b002d `i'b002e) if missingLS < 3
	drop missing*
}


foreach pheno in PA LS{
	rename KL_`pheno' `pheno'_2006
	rename ll_`pheno' `pheno'_2008
	rename ml_`pheno' `pheno'_2010
	rename nl_`pheno' `pheno'_2012
	rename ol_`pheno' `pheno'_2014
	rename pl_`pheno' `pheno'_2016
}


********************************************************
**************** Select variables, save ****************
********************************************************

keep hhid pn PA* LS* sex_* mob_* yob_* moi_* yoi_*
saveold "tmp/SWB_all.dta", replace version(12)
