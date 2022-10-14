clear all

local logfile="`1'"
local pheno_data="`2'"

log using `logfile', replace
display "$S_DATE $S_TIME"

set more off

import delimited `pheno_data'

**********************************************************
************************* ASTHMA *************************
**********************************************************
* https://www.opencodelists.org/codelist/opensafely/asthma-diagnosis/2020-04-15/

gen ASTHMA=0
* Got ctv3 codes from here https://www.opencodelists.org/codelist/opensafely/asthma-diagnosis/2020-04-15/download.csv , matched those two read_2 codes using all_lkps_maps_v3.xlsx downloaded from UKB data showcase (in PGI repo dropbox folder now, under 3.Phenotypes/UKB), took the union of ctv3 and read_2. Then examined the read_2 codes that were different than their ctv3 counterparts from the initial list to make sure they are not too broad (e.g. "did not attend asthma clinic" maps to just "did not attend" in read_2), dropped the ones that are. Also made sure meanings of newly added read_2 codes are not different in ctv3 (because I'm not separately going through the read_2 and read_3 columns below). Final list of codes below:
foreach i in "14B4." ".14B4" "173A." ".173A" "663N." ".663N" "663N0" "663N1" "663N2" "663O." ".663O" "663O0" "663P." ".663P" "663P0" "663P1" "663P2" "663Q." ".663Q" "663U." ".663U" "663V." ".663V" "663V0" "663V1" "663V2" "663V3" "663W." ".663W" "663d." ".663d" "663e." ".663e" "663e0" "663e1" "663f." ".663f" "679J." ".679J" "8791." ".8791" "8793." ".8793" "8794." ".8794" "8795." ".8795" "8796." ".8796" "8797." ".8797" "8798." ".8798" "8H2P." ".8H2P" "9OJ1." ".9OJ1" "9OJ3." ".9OJ3" "9OJ4." "9OJ5." ".9OJ5" "9OJ6." ".9OJ6" "9OJ7." ".9OJ7" "9OJ8." ".9OJ8" "9OJ9." ".9OJ9" "9OJA." ".9OJA" "9Q21." ".9Q21" "H3120" "H33.." "H330." "H3300" "H3301" "H330z" "H331." "H3310" "H3311" "H331z" "H332." "H33z." "H33z0" "H33z1" "H33z2" "H33zz" "H47y0" "Ua1AX" "X101t" "X101u" "X101x" "X101y" "X101z" "X1020" "X1021" "X1022" "X1023" "X1024" "X1025" "X1026" "X1027" "X1028" "X1029" "X102D" "X102G" "XE0YQ" "XE0YR" "XE0YS" "XE0YT" "XE0YU" "XE0YV" "XE0YW" "XE0YX" "XE0ZP" "XE0ZR" "XE0ZT" "XE2Nb" "XM0s2" "XM1Xb" "Xa0lZ" "Xa1hD" "Xa8Hn" "Xa9zf" "XaBAQ" "XaBU2" "XaBU3" "XaDvK" "XaDvL" "XaIIW" "XaIIX" "XaIIY" "XaIIZ" "XaINZ" "XaINa" "XaINb" "XaINc" "XaINd" "XaINf" "XaINg" "XaINh" "XaIOV" "XaIQ4" "XaIQD" "XaIQE" "XaIR3" "XaIRN" "XaIeq" "XaIer" "XaIfK" "XaIoE" "XaIu5" "XaIu6" "XaIuG" "XaIww" "XaJFG" "XaJYe" "XaKdk" "XaLIm" "XaLIn" "XaLIr" "XaLJS" "XaLJT" "XaLJU" "XaLPE" "XaNKw" "XaObi" "XaObj" "XaObk" "XaObl" "XaObm" "XaQHq" "XaQig" "XaQih" "XaQij" "XaR8K" "XaRFi" "XaRFj" "XaRFk" "XaRFl" "XaX3n" "XaXZm" "XaXZp" "XaXZs" "XaXZu" "XaXZx" "XaXa0" "XaY2V" "XaYZB" "XaYZh" "XaYb8" "XaYja" "XaYpF" "Xaa7B" "Xaa7Q" "Xabiu" "Xabj3" "Xac33" "Xac8r" "XacLz" "XacM0" "XacM1" "XacXj" "Xafdj" "Xafdy" "Xafdz" "Xaff0" "Y0017" "Y137e" "Y137f" "Y138a" "Y139b" "Y139c" "Y139d" "Y13a1" "Y13a2" "Y13a3" "Y13a4" "Y13a5" "Y13a7" "Y13a8" "Y13a9" "Y1b24" "Y1f9014Ok0" "173c." "173d." "178.." "1780." "1781." "1782." "1783." "1784." "1785." "1786." "1787." "1788." "1789." "178A." "178B." "1O2.." "388t." "388t0" "38DL." "38DT." "38DV." "661M1" "661N1" "663h." "663j." "663m." "663n." "663p." "663q." "663r." "663s." "663t." "663u." "663v." "663w." "663x." "663y." "66Y5." "66Y9." "66YA." "66YC." "66YE." "66YJ." "66YK." "66Yp." "66YP." "66Yq." "66YQ." "66Yr." "66YR." "66Ys." "66Yu." "66Yz5" "679J0" "679J1" "679J2" "8B3j." "8CMA0" "8CR0." "9NNX." "9OJ.." "9OJB." "9OJB0" "9OJB1" "9OJB2" "9OJC." "H302." "H333." "H334." "H335." "H3B.." {
        replace ASTHMA=2 if read_2=="`i'" | read_3=="`i'"
}

* Not doing control exclusions, gets too messy
bysort eid: egen ASTHMA_GP=max(ASTHMA)
**********************************************************

**********************************************************
************************** COPD **************************
**********************************************************
* H31.. (chronic bronchitis), H32.. (emphysema) except for senile/acute interstitial emphysema, obliterating fibrous bronchiolitis. Also searched for COPD, COAD, COLD, chronic obstructive pulmonary, chronic obstructive airway, chronic obstructove lung, chronic bronch, smoker's cough and included all relevant codes after inspection.
gen COPD=0
foreach i in ".14B3" ".66YB" ".66YD" ".66Yd" ".66Ye" ".66Yf" ".66Yg" ".66Yh" ".66YI" ".66Yi" ".66YL" ".66YM" ".66YS" ".66YT" ".679V" ".68C2" ".8CE6" ".8CR1" ".8H2R" ".9h5." ".9h51" ".9h52" ".9kf." ".9kf0" ".9kf1" ".9kf2" ".9N4W" ".9Oi." ".9Oi0" ".9Oi1" ".9Oi2" ".9Oi3" ".9Oi4" ".H37." ".H4.." ".H41." ".H411" ".H412" ".H413" ".H414" ".H415" ".H41Z" ".H42." ".H46." ".H47." ".H48." ".H4Z." ".H6Z2" "14B3." "14OX." "2126F" "38Dg." "661M3" "661N3" "66YB." "66YB0" "66YB1" "66YB2" "66YD." "66Yd." "66Ye." "66Yf." "66Yg." "66Yh." "66YI." "66Yi." "66YL." "66YM." "66YS." "66YT." "66Yz1" "66Yz2" "679V." "8BMa0" "8BMW." "8CE6." "8CeD." "8CMV." "8CMW5" "8CR1." "8H2R." "8Hkw." "8I610" "8IEy." "8IEZ." "9e03." "9h5.." "9h51." "9h52." "9kf.." "9kf0." "9kf1." "9kf2." "9N4W." "9NgP." "9Nk70" "9Oi.." "9Oi0." "9Oi1." "9Oi2." "9Oi3." "9Oi4." "H31.." "H310." "H3100" "H3101" "H310z" "H311." "H3110" "H3111" "H311z" "H312." "H3120" "H3121" "H3122" "H3123" "H312z" "H313." "H31y." "H31y0" "H31y1" "H31yz" "H31z." "H32.." "H320." "H3200" "H3201" "H3202" "H3203" "H320z" "H321." "H322." "H32y0" "H32y2" "H32yz" "H32z." "H36.." "H37.." "H38.." "H39.." "H3A.." "H3B.." "H3y.." "H3y0." "H3y1." "H3z.." "Hyu30" "Hyu31" "X101i" "X101j" "X101l" "X101n" "X101o" "Xa35l" "Xaam4" "XaavA" "Xabwy" "Xac33" "Xac8s" "XaEIV" "XaEIW" "XaEIY" "Xafex" "Xafhe" "XafhZ" "Xafir" "Xafit" "Xafiu" "XaIes" "XaIet" "XaIND" "XaIql" "XaIQT" "XaIRO" "XaIu7" "XaIu8" "XaIUt" "XaJ4k" "XaJ4l" "XaJ4R" "XaJDW" "XaJFu" "XaJlS" "XaJlT" "XaJlU" "XaJlV" "XaJlW" "XaJlY" "XaJYf" "XaK8Q" "XaK8R" "XaK8S" "XaK8U" "XaKv8" "XaKv9" "XaKzy" "XaLJz" "XaLqj" "XaN4a" "XaPio" "XaPls" "XaPlu" "XaPZH" "XaPzu" "XaRCG" "XaRCH" "XaRFy" "XaW9D" "XaX3c" "XaXCa" "XaXCb" "XaXnt" "XaXzy" "XaY05" "XaY0w" "XaYbA" "XaYuZ" "XaYZO" "XaZ6U" "XaZ8t" "XaZd1" "XaZee" "XaZoz" "XaZp7" "XaZp8" "XaZp9" "XE0tj" "XE0YM" "XE0YN" "XE0YP" "XE0ZL" "XE0ZN" "XE2Pp" {
        replace COPD=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls: H3... Chronic obstructive pulmonary disease (too broad, includes asthma and other things), H32... Emphysema and H32y Other emphysema (too broad, includes senile/acute interstitial emphysema)
foreach i in "H3..." "H32.." "H32y."  {
        replace COPD=1 if COPD==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen COPD_GP=max(COPD)
**********************************************************


**********************************************************
************************* MIGRAINE ***********************
**********************************************************
gen MIGRAINE=0
foreach i in "1474" "14740" ".F351" "F260." "F261." "F2610" "F261z" "F26y0" "X007J" "X007K" "X007L" "X007M" "X007N" "XaXkv" "0.1474" ".8B6N" ".F35." "8B6N." "F2623" "0.1967" "1967" "F2622" "F2624" "F262z" "F26y." "F26y0" "F26y1" "F26y2" "F26y3" "F26yz" "F26z." "Fyu53" "R090D" "X007O" "X007S" "XaJLO" "K584." "Xa07H" "XaXkr" "F2628"{
        replace MIGRAINE=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls: sick headache, neuralgic migraine/migrainous neuralgia/cluster headache (because they're coded the same), unspecified migraine and migraine variants (because they include cluster headache, tension headaches, etc which are not classified as migraine in ICD9/10)
foreach i in "F2611" "F2621" ".F35Z" "F26.." "F262." "F2620" "F2625" "XE187" {
        replace MIGRAINE=1 if MIGRAINE==0 & (read_2=="`i'" |  read_3=="`i'")
}
bysort eid: egen MIGRAINE_GP=max(MIGRAINE)
**********************************************************

**********************************************************
************************ ECZEMA **************************
**********************************************************
gen ECZEMA=0
foreach i in ".L23." ".L232" ".L233" ".L234" ".L23Z" "M111." "M112." "M113." "M114." "M115." "M116." "M117." "M11z." "X505e" "X505f" "X505M" "X505N" "X505O" "X505P" "X505Q" "X505R" "X505S" "X505T" "X505U" "X505V" "X505Z" "X506c" "XaBsL" "XE1C6" {
        replace ECZEMA=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in ".14F1" ".26C4" ".F5C4" ".L2.." ".L27." ".L28." ".L2Z." ".L318" ".L342" ".L353" "14F1." "26C4." "38Vp." "8HTu." "F4D30" "F5024" "M070." "M07y." "M07z." "M1..." "M102." "M11.." "M110." "M11A." "M12z0" "M12z2" "M12z3" "M12z4" "M1535" "M1536" "M173." "M1831" "M1832" "M21A." "M2522" "m5..." "Myu2." "Myu22" "X00iS" "X00iT" "X30Cp" "X40Fx" "X5040" "X504r" "X505K" "X505L" "X5061" "X506d" "X506O" "X506P" "X506Q" "X506T" "X506U" "X50AK" "X78z3" "Xa0p8" "Xa1dl" "Xa4jb" "Xa7lb" "Xa7lZ" "XaBml" "XaBmm" "Xaepv" "XaINK" "XaINM" "XaL2Q" "XaQfn" "XaY4a" "XaY4o" "XaY4Z" "XE1Av" "XE1CI" "XE1CW" "XM1GA" "XM1Pa"  {
        replace ECZEMA=1 if ECZEMA==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen ECZEMA_GP=max(ECZEMA)
**********************************************************

**********************************************************
*********************** HAYFEVER *************************
**********************************************************
gen HAYFEVER=0
foreach i in "14B1." "H17.." "H170." "H171." "H1710" "H1711" "H172." "H17z." "Hyu20" "Hyu21" "14M4." ".14B1" ".14M4" ".H28." ".H290" ".PD.." "X00l3" "X00l4" "X00l5" "X00l6" "X00l7" "X00l8" "X00l9" "X00lA" "X00lB" "Xa0lX" "Xa7lL" "Xa7lM" "XaIO5" "XaIpW" "XaIRZ" "XaOb5" "XaOj7" "XE0Y5" "XE0Y6" "XE0Y7" "XE2QI" {
        replace HAYFEVER=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in "H1..." "Hyu2." "H330." "H3300" ".H29." "X1020" "X70vw" "Xa0Li" "Xa0Lj" "XE0Z5" "XM1Md" {
        replace HAYFEVER=1 if HAYFEVER==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen HAYFEVER_GP=max(HAYFEVER)
**********************************************************

**********************************************************
*********** INFLAMMATORY BOWEL DISEASE *******************
**********************************************************
gen IBD=0
foreach i in ".I51." ".I511" ".I512" ".I51Z" ".I52." "8CA4W" "8Cc5." "idD.." "J08z9" "J40.." "J400." "J4000" "J4001" "J4002" "J4003" "J4004" "J4005" "J400z" "J401." "J4010" "J4011" "J4012" "J401z" "J402." "J40z." "J41.." "J410." "J4100" "J4101" "J4102" "J4103" "J4104" "J410z" "J411." "J412." "J413." "J41y." "J41y0" "J41yz" "J41z." "J436." "J4360" "J4361" "J438." "J4z6." "Jyu40" "Jyu41" "N0310" "N0311" "N0453" "N0454" "X20Pq" "X300J" "X301b" "X302r" "X302t" "X302u" "X302y" "X3030" "X303k" "X303x" "X303y" "X3050" "X7021" "X702C" "Xa0lh" "XaaSE" "XabAl" "XaK6C" "XaK6D" "XaK6E" "XaXla" "XaYzX" "XaZ2j" "XE0ae" "XE0af" "XE0ag" "XE0cZ" "XE2QL" "XM0bn" "XM1RP" {
        replace IBD=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in ".14C4" ".I5.." ".I6ZZ" "J4..." "J41y1" "J43.." "J437." "J43z." "J4z.." "J4z1." "J4z2." "J4z3." "J4z4." "J4z5." "J4zz." "J574C" "J6A.." "Jyu4." "Jyu42" "X302j" "X302q" "X302s" "X303e" "X303f" "X303j" "X3047" "X304s" "X30BN" "X70Ft" "Xa7nY" "XaB5x" "XaB5z" "XaC1C" "XaK6F" "XC01B" "XC0tD" "XE0aj" "XE0an" "XE0ao" "XE0cX" "XE0d7" "XE0qC" {
        replace IBD=1 if IBD==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen IBD_GP=max(IBD)
**********************************************************


**********************************************************
******************** TYPE-II DIABETES ********************
**********************************************************
gen T2D=0
foreach i in ".66Ao" ".C21." ".C2A." ".C2D." ".C2D0" "66Ao." "66At1" "C1001" "C1011" "C1021" "C1031" "C1041" "C1051" "C1061" "C1071" "C1072" "C1074" "C109." "C1090" "C1091" "C1092" "C1093" "C1094" "C1095" "C1096" "C1097" "C1099" "C109A" "C109B" "C109C" "C109D" "C109E" "C109F" "C109G" "C109H" "C109J" "C109K" "C10F." "C10F0" "C10F1" "C10F2" "C10F3" "C10F4" "C10F5" "C10F6" "C10F7" "C10F8" "C10F9" "C10FA" "C10FB" "C10FC" "C10FD" "C10FE" "C10FF" "C10FG" "C10FH" "C10FJ" "C10FK" "C10FL" "C10FM" "C10FN" "C10FP" "C10FQ" "C10FR" "C10P1" "C10y1" "C10z1" "L1806" "L180B" "X40J5" "X40J6" "Xa2hA" "Xaagf" "XaELQ" "XaEnp" "XaEnq" "XaF05" "XaFmA" "XaFn7" "XaFn8" "XaFn9" "XaFWI" "XaIfG" "XaIfI" "XaIrf" "XaIzQ" "XaIzR" "XaJQp" "XaKyX" "XaMhK" "XaX3q" "XaXZR" "XE10F" "XE12A" "XM19j" {
        replace T2D=2 if read_2=="`i'" | read_3=="`i'"
}

* Not doing control exclusions, too many unspecified diabetes categories
bysort eid: egen T2D_GP=max(T2D)
**********************************************************

**********************************************************
********************** ALLERGY - CAT *********************
**********************************************************
gen ALLERGYCAT=0
foreach i in "X00l3" "XaIO5" "x001l" "x00CK" "x00CT" "x00Cc" "x00Cl" ".14M4" "14M4."{
        replace ALLERGYCAT=2 if read_2=="`i'" | read_3=="`i'"
}

* Exclude from controls: dander allergy, animal hair, animal, allergic rhinitis due to other allergens
foreach i in "H171." "Xa7lM" "XaIpW" ".H290" "H1710"{
        replace ALLERGYCAT=1 if ALLERGYCAT==0 & (read_2=="`i'" | read_3=="`i'")
}

bysort eid: egen ALLERGYCAT_GP=max(ALLERGYCAT)

**********************************************************

**********************************************************
******************* ALLERGY - POLLEN *********************
**********************************************************
gen ALLERGYPOLLEN=0
foreach i in "H170." "H330." "X1020" "XE0ZP" "XE2QI" "Xa7lL" "XaObk" "XaOj7" "c91.." "c911." "c912." "c913." "c914." "c915." "c916." "c917." "c918." "c919." "c91a." "c91b." "c91c." "c91d." "c91e." "c91f." "c91g." ".H431" "1781." "Xa0lX" "x0563" "x0564" "x0565" "x0566" ".H28." "14M2." ".14M2" "x00CR" "x00Ca" "x00Cj" "x00Ct" "x00CM" "x00CQ" "x00CV" "x00CZ" "x00Ce" "x00Ci" "x00Co" "x00Cs" ".14B1" "14B1." "Hyu20" "1781." "X71bE" {
        replace ALLERGYPOLLEN=2 if read_2=="`i'" | read_3=="`i'"
}

bysort eid: egen ALLERGYPOLLEN_GP=max(ALLERGYPOLLEN)

**********************************************************

**********************************************************
********************* ALLERGY - DUST *********************
**********************************************************
gen ALLERGYDUST=0
foreach i in "X00l6" "X00l7" "XaIRZ" "c92.." "c921." "c922." "c923." "c924." "c925." "c926." "c927." "c928." "x005O" "x00CO" "x00CX" "x00Cg" "x00Cq" ".PD.." "x005P" "x00CP" "x00CY" "x00Ch" "x00Cr" "X71bD"{
        replace ALLERGYDUST=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude from controls: allergic rhinitis due to other allergens
foreach i in "H171." {
        replace ALLERGYDUST=1 if ALLERGYDUST==0 & (read_2=="`i'" | read_3=="`i'")
}
bysort eid: egen ALLERGYDUST_GP=max(ALLERGYDUST)

**********************************************************


**********************************************************
********************* NEARSIGHTEDNESS ********************
**********************************************************
gen NEARSIGHTED=0
foreach i in "F4021" "F471." "X00cl" "X00g3" "X00g4" "X00g5" "X00g6" "X75oH" "X75oI" "X75oJ" "XaE5J" ".F562" "F4025" "F4710"{
        replace NEARSIGHTED=2 if read_2=="`i'" | read_3=="`i'"
}
bysort eid: egen NEARSIGHTED_GP=max(NEARSIGHTED)
**********************************************************


**********************************************************
******************** ANOREXIA NERVOSA ********************
**********************************************************
gen ANOREX=0
foreach i in "1467." "E271." "X00Sz" ".1467" ".E49." "Eu500" "Eu501" {
        replace ANOREX=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in "1612." "R030." "R030z" "X76cG" "XE24f" "XM07X" "XM0aM" "XM0aN" ".1612" ".R30." ".R30Z" "E275z" "Eu50." "Eu50y" "Eu50z" "X00Sx" "XE1Zq" "XE1bQ" "Xa2hW" "XaEC2" "XaIwi" "XaPBE" "XaWRh" "XaX4J" "XaX4K" "XaX4L" "XaX4M" "XafF2" "XafjT" ".8HTN" ".9Nk9" ".E4D4" "38Do." "38Do0" "38Do1" "38Do2" "38Do3" "8HTN." "9Nk9." {
        replace ANOREX=1 if ANOREX==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen ANOREX_GP=max(ANOREX)
**********************************************************



**********************************************************
************************ BIPOLAR *************************
**********************************************************
gen BIPOLAR=0
foreach i in ".146D" ".212V" ".E22." "146D." "212V." "E111." "E1110" "E1111" "E1112" "E1113" "E1114" "E1115" "E1116" "E111z" "E114." "E1140" "E1141" "E1142" "E1143" "E1144" "E1145" "E1146" "E114z" "E115." "E1150" "E1151" "E1152" "E1153" "E1154" "E1155" "E1156" "E115z" "E116." "E1160" "E1161" "E1162" "E1163" "E1164" "E1165" "E1166" "E116z" "E117." "E1170" "E1171" "E1172" "E1173" "E1174" "E1175" "E1176" "E117z" "E11y." "E11y0" "E11y1" "E11y3" "E11yz" "Eu31." "Eu310" "Eu311" "Eu312" "Eu313" "Eu314" "Eu315" "Eu316" "Eu317" "Eu318" "Eu319" "Eu31y" "Eu31z" "X00SM" "X00SN" "XaCHo" "XaIx7" "XaMwc" "XaY1Y" "XE1aQ" "XE1ZX" "ZV111" "XaB95" "Eu332" "Eu333" {
        replace BIPOLAR=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in ".13Y3" ".212T" ".E22." ".E22Z" "13Y3." "212T." "E11.." "Eu30." "Eu332" "Eu333" "XaA6j" "XafEv" "XaLIa" "XE1aQ" {
        replace BIPOLAR=1 if BIPOLAR==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen BIPOLAR_GP=max(BIPOLAR)
**********************************************************



**********************************************************
******************** SCHIZOPHRENIA ***********************
**********************************************************
gen SCZ=0
foreach i in ".1464" ".212W" ".E21." ".E211" ".E212" ".E213" ".E214" ".E21Z" "1464." "212W." "E10.." "E100." "E1000" "E1001" "E1002" "E1003" "E1004" "E1005" "E100z" "E101." "E1010" "E1011" "E1012" "E1013" "E1014" "E1015" "E101z" "E102." "E1020" "E1021" "E1022" "E1023" "E1024" "E1025" "E102z" "E103." "E1030" "E1031" "E1032" "E1033" "E1034" "E1035" "E103z" "E104." "E106." "E107." "E1070" "E1071" "E1072" "E1073" "E1074" "E1075" "E107z" "E10y." "E10y0" "E10y1" "E10yz" "E10z." "Eu20." "Eu200" "Eu201" "Eu202" "Eu203" "Eu204" "Eu205" "Eu206" "Eu20y" "Eu20z" "Eu25." "Eu250" "Eu251" "Eu252" "Eu25y" "Eu25z" "X00S8" "X00SF" "X00SG" "X00SH" "X00SI" "XaMwd" "XE1aM" "XE1aO" "XE1Xx" "XE1ZM" "XE2b8" "XE2un" "XE2uT" "ZV110" {
        replace SCZ=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls
foreach i in "13Y2." "E14z." "Eu2.." "X00S7" "X00TP" "Xa1c7" "XaLIa" ".13Y2" ".212T" ".E2Z2" "212T." "Eu845" "X761M" {
        replace SCZ=1 if SCZ==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen SCZ_GP=max(SCZ)
**********************************************************


**********************************************************
*************************** ADHD *************************
**********************************************************
gen ADHD=0
foreach i in ".E4D6" ".E4D7" ".E4H." ".E4H1" "6A61." "8BPT." "8BPT0" "8BPT1" "9Ngp." "9Ngp0" "9Ngp1" "9Ol8." "9Ol9." "9OlA." "E2E.." "E2E0." "E2E00" "E2E01" "E2E0z" "E2E1." "E2E2." "E2Ey." "E2Ez." "Eu90." "Eu900" "Eu901" "Eu902" "Eu90y" "Eu90z" "Eu9y7" "Ub1Tu" "Ub1Tv" "Ub1Tw" "Ub1Tx" "Ub1Ty" "Ub1Tz" "Xa1bR" "Xaae7" "Xaae8" "Xaae9" "XaaZA" "XaaZj" "XaaZk" "XabmL" "XaQxu" "XaQzA" "XaQzB" "XaVw9" "XE2Q6" "ZV405" {
        replace ADHD=2 if read_2=="`i'" | read_3=="`i'"
}
* Exclude these from controls: Childhood conduct disorder, Hyperactive behaviour
foreach i in "X00TU" "X764X" ".1P00" "1P00." {
        replace ADHD=1 if ADHD==0 & (read_2=="`i'" |  read_3=="`i'")
}

bysort eid: egen ADHD_GP=max(ADHD)
**********************************************************


bysort eid: gen i=_n
drop if i>1
keep eid *_GP
ren eid n_eid
save "tmp/gp_clinical.dta", replace

foreach pheno in *GP {
        tab `pheno'
}

log close