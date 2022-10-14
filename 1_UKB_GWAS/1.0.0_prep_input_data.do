clear all

local crosswalk="/disk/genetics2/ukb/orig/UKBv3/crosswalk/ukb_imp_chr1_22_v3.ukb11425_imp_chr1_22_v3_s487395.crosswalk"
local partition_data="/disk/genetics/PGS/PGI_Repo/derived_data/1_UKB_GWAS/tmp/IDs_assignPartition_ordered.txt"
local pheno_data_1="/disk/genetics2/ukb/orig/UKBv2/pheno/9690_update/ukb9690_update_jun14_2018.dta"
local pheno_data_2="/disk/genetics/ukb/data/41490_alc_cannabis_menarche_myopia/ukb41490.dta"
local pheno_data_3="/disk/genetics2/ukb/orig/UKBv2/pheno/52008/ukb52008.dta"
local covar_data="/disk/genetics2/ukb/orig/UKBv2/linking/ukb_sqc_v2_combined_header.dta"
local logfile="/disk/genetics/PGS/PGI_Repo/code/1_UKB_GWAS/1.0.0_prep_input_data.do.log" 
local withdrawn="/disk/genetics/ukb/data/withdrawals/ukb_withdrawn_ind_20220222.csv"

local crosswalk="`1'"
local partition_data="`2'"
local pheno_data_1="`3'"
local pheno_data_2="`4'"
local pheno_data_3="`5'"
local covar_data="`6'"
local withdrawn="`7'"
local logfile="`8'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off
set maxvar 32000

use `pheno_data_1'
merge 1:1 n_eid using `pheno_data_2', nogen
keep n_eid      /// ID
    n_22010_*   /// Bad genotype
    n_34_0_0    /// Birth year
    n_20262_0_0 /// Myopia diagnosis (from pheno_data_2)

merge 1:1 n_eid using `pheno_data_3', nogen

gen IID=n_eid
merge 1:1 IID using `covar_data', nogen keep(match)

keep n_eid 			/// ID
    n_22010_*   /// Bad genotype
    n_34_0_0    /// Birth year
    PC*             /// PCs (from covar)
    Sex             /// Sex (from covar)
    InferredGender  /// Genetic sex (from covar)
    SubmittedGender /// Reported sex (from covar)
    Batch           /// Genotyping batch (from covar)
    hetmissingoutliers /// (from covar)
    putativesexchromosomeaneuploidy /// (from covar)
	n_21000_* 		/// Ethnicity 				/ 0-2
    n_21003_* 		/// age at assesment visit	/0-3
    n_22200_0_0 	/// birth year
    n_22027_0_0 	/// bad genotype?
    n_22001_0_0 	/// Genetic sex
	s_41270_* 		/// ICD10 summary
	s_41271_* 		/// ICD9 summary
	s_41202_* 		/// ICD10 primary			/ 0 0-78
	s_41204_* 		/// ICD10 secondary			/ 0 0-187
	s_41203_* 		/// ICD9 primary			/ 0-27
	s_41205_* 		/// ICD9 secondary			/ 0-29
	s_40001_* 		/// ICD10 death primary		/ 0-1
	s_40002_* 		/// ICD10 death secondary	/ 0 1-14, 1 1-9
	s_40006_* 		/// ICD10 type of cancer	/ 0-21
    s_40013_* 		/// ICD9 type of cancer 	/ 0-14
	s_41272_* 		/// Operations OPCS4 		/ 0-123
	s_41273_* 		/// Operations OPCS3 		/ 0-15
	n_20001_* 		/// Self-reported cancer 	/ 0 0-5, 1 0-3, 2 0-4, 3 0-3
	n_20002_* 		/// Self-reported non-cancer illness  /0 0-28, 1 0-15, 2 0-33, 3 0-18
	n_20003_*		/// Self-reported medications
	n_20004_*		/// Self-report operations	/ 0 0-31, 1 0-14, 2 0-17, 3 0-9
	n_22126_* 		/// doctor diagnosed hayfever or allergic rhinitis
	n_6150* 		/// Other diagnoses (heart problems) 0 0-3, 1 0-3, 2 0-3
	n_6152* 		/// Other diagnoses (Hayfever)  0 0-4, 1 0-3, 2 0-4, 3 0-2
	n_6153_* 		/// Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones	0 0-3, 1 0-3, 2 0-3, 3 0-2
	n_6177_* 		/// Medication for cholesterol, blood pressure or diabetes	0-3 x 0-2
	n_2453_* 		/// ever doctor diagnoed with cancer 3
    n_30780_*       /// LDL cholesterol
    ///n_30781_*       /// LDL cholesterol assay date
    n_30782_*       /// LDL cholesterol aliquot
    n_30760_*       /// HDL cholesterol
    ///n_30761_*       /// HDL cholesterol assay date
    n_30762_*       /// HDL cholesterol aliquot
    n_30690_*       /// Total cholesterol
    ///n_30691_*       /// Total cholesterol assay date
    n_30692_*       /// Total cholesterol aliquot
    n_30870_*       /// Tryglycerids
    ///n_30871_*       /// Tryglycerids assay date
    n_30872_*       /// Tryglycerids aliquot
    ts_3166_*       /// Time blood sample collected
    n_74_*          /// Fasting time
    n_30897_*       /// Estimated sample dilution factor
    n_54_*          /// Assessment center
    n_55_*          /// Month of assessment
	n_2178* 		/// Self-rated Health 3
	n_2000_* 		/// Worry too long after embarrassment 0-3
    n_2010_* 		/// Suffer from 'nerves'	 0-3
    n_2020_* 		/// Loneliness, isolation 0-3
    n_2030_* 		/// Guilty feelings	 0-3
	n_2040_* 		/// Risk 0-3
	n_2050_* 		/// Frequency of depressed mood in last 2 weeks	
	n_2060_* 		/// Frequency of unenthusiasm
	n_2070_* 		/// Frequency of tenseness / restlessness in last 2 weeks
	n_2080_* 		/// Depressive symptoms 3
	n_4526*  		/// Happiness / Subjective Wellbeing 3 
	n_4537_* 		/// Work satisfaction  0-3
	n_4559_* 		/// Family satisfaction 0-3
    n_4570_* 		/// Friendship satisfaction 0-3
    n_4581_* 		/// Financial satisfaction 0-3
	n_20122_0_0 	/// Bipolar status
	n_20126_0_0 	/// Bipolar and MDD status
	n_20127_0_0 	/// Neuroticism score 
	n_20544_* 		/// Mental health problems ever diagnosed by a professional 1-16
	n_130875_0_0	/// Source of report of F20 (schizophrenia)	
	n_130885_0_0	/// Source of report of F25 (schizoaffective disorders)
	n_130893_0_0 	/// Source of bipolar affective disorder
	n_20499_0_0		/// Ever sought or received professional help for mental distress
	n_20500_0_0		/// Ever suffered mental distress preventing usual activities
	n_4079_* 		/// Diastolic blood pressure, automated reading 0-3 x 0-1
    n_94_* 			/// Diastolic blood pressure, manual reading 0-3 x 0-1
    n_4080_* 		/// Systolic blood pressure, automated reading 0-3 x 0-1
    n_93_* 			/// Systolic blood pressure, manual reading 0-3 x 0-1
	n_20117_* 		/// EVERDRINK status (derived)
    n_20403_0_0  	/// Amount of alcohol drunk on a typical drinking day 1
    n_20405_0_0  	/// Ever had known person concerned about, or recommend reduction of, alcohol consumption 1
    n_20407_0_0  	/// Frequency of failure to fulfil normal expectations due to drinking alcohol in last year 1
    n_20408_0_0  	/// Frequency of memory loss due to drinking alcohol in last year 1
    n_20409_0_0  	/// Frequency of feeling guilt or remorse after drinking alcohol in last year 1
    n_20411_0_0  	/// Ever been injured or injured someone else through drinking alcohol 1
    n_20412_0_0  	/// Frequency of needing morning drink of alcohol after heavy drinking session in last year 1
    n_20413_0_0  	/// Frequency of inability to cease drinking in last year 1
    n_20414_0_0  	/// Frequency of drinking alcohol 1
    n_20416_0_0  	/// Frequency of consuming six or more units of alcohol 1
	n_1558_* 		/// drinking frequency 0-3
    n_1568_* 		/// red wine per week 0-3
    n_1578_* 		/// white wine and champagne 0-3
    n_1588_* 		/// beer and cider per week 0-3
    n_1598_* 		/// spirits per week 0-3
    n_1608_* 		/// fortified wine per week 0-3
    n_5364_* 		/// other alcoholic drinks per week LESS PEOPLE 0-3
    n_4407_* 		/// red wine per month 0-3
    n_4418_* 		/// white and champagne 0-3
    n_4429_* 		/// beer
    n_4440_* 		/// spirits 0-3
    n_4451_* 		/// fortified wine 0-3
    n_4462_* 		/// other 0-3
	n_1239_* 		/// current smoking 0-3
    n_1249_* 		/// past smoking 0-3
    n_2644_* 		/// light smoking 0-3
	n_2867_*		/// Age started smoking in former smokers
	n_2887_* 		/// Number of cigarettes previously smoked daily
	n_3436_*  		/// Age started smoking in current smokers
	n_3456_* 		/// Number of cigarettes currently smoked daily (current cigarette smokers)	
	n_6183_* 		/// Number of cigarettes previously smoked daily (current cigar/pipe smokers)	  0-3
	n_20160* 		/// Ever smoked 3
    n_20116_* 		/// smoking status 0-3
	n_20453_0_0 	/// cannabis
	n_2139_* 		/// Age at first sex 3
	n_2714_* 		/// Age at menarche
	n_2754* 		/// Age at first birth 3
    n_2734* 		/// Number of children (women) 3
    n_2405* 		/// Number of children (men) 3
    n_3872* 		/// Age at first birth single child 0-3
	n_50_* 			/// Height 3
	n_23104_*		/// BMI 3
    n_1180_* 		/// Morning person 3
	n_1200_* 		/// Insomnia 3
	n_6160* 		/// Religious group activity 0-3 x 0-4
    n_6138_* 		/// EA 0-3 x 0-5 
    n_845_* 		/// Age left schooling 0-2
    n_20016_* 		/// CP touchscreen 0-3
    n_20191_0_0 	/// CP web-based 1
    n_2207_* 		/// Wears glasses or contact lenses 0-3
	n_6147_* 		/// Reason for glasses/contact lenses 0-3
    n_20262_0_0     /// Myopia diagnosis 
    n_120016_0_0    /// Ever had migraine
    n_22126_0_0     /// Doctor diagnosed hayfever or allergic rhinitis
    n_22127_0_0     /// Doctor diagnosed asthma
    n_22128_0_0     /// Doctor diagnosed emphysema
    n_22129_0_0     /// Doctor diagnosed chronic bronchitis
    n_22130_0_0     /// Doctor diagnosed COPD
    n_130709_0_0    /// Source of report of E11 (non-insulin-dependent diabetes mellitus)
    n_131629_0_0    /// Source of report of K51 (ulcerative colitis)
    n_131627_0_0    /// Source of report of K50 (crohn's disease [regional enteritis])
    n_20107_*       /// Illnesses of father
    n_20110_*       /// Illnesses of mother
    n_20111_*       /// Illnesses of siblings
    n_42020_*       ///	Date of alzheimer's disease report	
    n_42021_*       ///	Source of alzheimer's disease report
    n_2946_*        /// Father's age
    n_1807_*        /// Father's age at death
    n_1845_*        /// Mother's age
    n_3526_*        /// Mother's age at death


    
*** QC ***
* Missing sex
drop if Sex == 0 
drop if missing(Sex)

* Sex mismatch
drop if InferredGender != SubmittedGender

* Missing batch
drop if missing(Batch)

* Bad genotypes, second release
drop if hetmissingoutliers == 1
drop if missing(PC1-PC40)
drop if n_eid < 0 
drop if putativesexchromosomeaneuploidy == 1

* Bad genotypes, initial release
drop if n_22010_0_0 == 1 

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

gen BYEAR2 = BYEAR*BYEAR
gen BYEAR3 = BYEAR2*BYEAR
gen SexxBYEAR = Sex*BYEAR
gen SexxBYEAR2 = Sex*BYEAR2
gen SexxBYEAR3 = Sex*BYEAR3

save "tmp/pgi_repo.dta", replace


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
save "tmp/pgi_repo.dta", replace
