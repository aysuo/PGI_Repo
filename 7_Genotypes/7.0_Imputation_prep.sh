#!/bin/bash

source $PGI_Repo/code/paths
pathGf=$1
refpanel=$2

if ! [[ -s $PGI_Repo/code/7_Genotypes/7.0.0_HRC-1000G-check-bim.pl ]]; then
    wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip -O $PGI_Repo/code/7_Genotypes/HRC-1000G-check-bim-v4.3.0.zip
    unzip $PGI_Repo/code/7_Genotypes/HRC-1000G-check-bim-v4.3.0.zip -d $PGI_Repo/code/7_Genotypes/
    mv $PGI_Repo/code/7_Genotypes/HRC-1000G-check-bim.pl $PGI_Repo/code/7_Genotypes/7.0.0_HRC-1000G-check-bim.pl
    rm $PGI_Repo/code/7_Genotypes/HRC-1000G-check-bim-v4.3.0.zip $PGI_Repo/code/7_Genotypes/LICENSE.txt
fi

if [[ $refpanel == "HRC" ]]; then
    if ! [[ -s $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab ]]; then
        wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz -O $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
        gunzip $PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
    fi
    refdata=$PGI_Repo/original_data/ref_data/HRC_imputation_qc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
    
elif [[ $refpanel == "TOPmed" ]]; then
    if ! [[ -s $PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/PASS.Variants.TOPMed_freeze3a_hg19_dbSNP.tab ]]; then
        curl 'https://bravo.sph.umich.edu/freeze3a/hg19/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _ga=GA1.2.1019045455.1673007919; _gid=GA1.2.594407710.1673007919; remember_token="aysuokbay@gmail.com|fa2416580a6d934deea0a8c00744d815816f223bd1b34823a38f4523ae4dfc1d23676bab50668a2ce6bdb7ee8764af7da61f24c6f1b4196f74f1a74dcc603133"; _gat_gtag_UA_73910830_2=1; _gat=1' --compressed > $PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/ALL.TOPMed_freeze3a_hg19_dbSNP.vcf.gz
        wget https://www.well.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip -O $PGI_Repo/code/7_Genotypes/CreateTOPMed.zip
        unzip $PGI_Repo/code/7_Genotypes/CreateTOPMed.zip -d $PGI_Repo/code/7_Genotypes/
        mv $PGI_Repo/code/7_Genotypes/CreateTOPMed.pl $PGI_Repo/code/7_Genotypes/7.0.1_CreateTOPMed.pl
        $PGI_Repo/code/7_Genotypes/7.0.1_CreateTOPMed.pl -i $PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/ALL.TOPMed_freeze3a_hg19_dbSNP.vcf.gz -o $PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/PASS.Variants.TOPMed_freeze3a_hg19_dbSNP.tab.gz
        gunzip $PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/PASS.Variants.TOPMed_freeze3a_hg19_dbSNP.tab.gz
        rm $PGI_Repo/code/7_Genotypes/CreateTOPMed.zip $PGI_Repo/code/7_Genotypes/LICENSE.txt
    fi
    refdata=$PGI_Repo/original_data/ref_data/TOPmed_imputation_qc/PASS.Variants.TOPMed_freeze3a_hg19_dbSNP.tab
fi


getFreq(){
    plink2 --bfile ${pathGf} \
    --freq \
    --out ${pathGf}_frq
}

HRC_1000G_check_bim(){
    perl  $PGI_Repo/code/7_Genotypes/7.0.0_HRC-1000G-check-bim.pl \
        -b ${pathGf}.bim \
        -f ${pathGf}_frq.afreq \
        -h \
        -r ${refdata}
    
    # "plink" calls plink1.7 on the server, change that to plink1.9
    sed -i 's/plink/plink1.9/g' Run-plink.sh
    sh Run-plink.sh | tee Run-plink.sh.log
}

convert2vcf(){
    for chr in {1..22}; do
        plink1.9 --bfile ${pathGf}-updated-chr$chr --recode vcf --out ${pathGf}_qc_chr$chr
    done
}

sortVcf(){
    for chr in {1..22}; do
        $vcfsort ${pathGf}_qc_chr$chr.vcf | bgzip -c > ${pathGf}_qc_chr${chr}_sorted.vcf.gz
    done
}

checkVcf(){
    for chr in {1..22}; do
	    $python $checkvcf/checkVCF.py \
            -r $checkvcf/hs37d5.fa \
            -o ${pathGf}_qc_chr${chr}_vcfcheck ${pathGf}_qc_chr${chr}_sorted.vcf.gz
    done
}


################################################################################################

HRCqc(){
    dirGf=$(echo $pathGf | rev | cut -d"/" -f1 --complement | rev)

    cd $dirGf
    getFreq 
    HRC_1000G_check_bim 
    convert2vcf 
    sortVcf 
    checkVcf
}


HRCqc
