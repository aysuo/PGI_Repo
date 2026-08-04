# PGI_Repo
This repository contains the complete pipeline for the creation of polygenic indexes (PGIs) for the second release of the SSGAC Polygenic Index Repository, an initiative that makes PGIs for a wide range of traits available for a number of datasets. The Repository currently contains PGIs for [61 phenotypes](https://www.dropbox.com/scl/fi/z8oen27v7n74y3fqv9ce7/Datasets-participating-in-the-Repository.pdf?rlkey=o6abg8ik3e1ezmrivp2i7y3ie&st=1hthmctb&dl=0) in [22 datasets](https://www.dropbox.com/scl/fi/z8oen27v7n74y3fqv9ce7/Datasets-participating-in-the-Repository.pdf?rlkey=o6abg8ik3e1ezmrivp2i7y3ie&st=1hthmctb&dl=0). Please see [Alemu et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.05.14.653986v2)  and the [User Guide](https://www.dropbox.com/scl/fo/o1l5c12iyztay4evlxn0l/AAl-n6cDJl0xIMHDmM2VRSk?rlkey=texlne2to5gqbqxweovzwozef&st=biahlpsl&dl=0) for a detailed description of the pipeline. 

The Repository will be updated regularly with additional PGIs and datasets. If you are interested in participating in the Repository, please reach out to contact@ssgac.org. 

## Frequently Asked Questions (FAQs)
For a less technical description of the paper and of how PGIs should—and should not—be interpreted and used, see these [frequently asked questions](https://www.dropbox.com/scl/fi/jvof57tbyhp33g34cuzic/FAQ-Polygenic-Index-Repository.pdf?rlkey=9xc3cxeqjtd9iyqkcue11uuap&st=rhexfjy4&dl=0).

## PGI Access Procedures
PGIs in the participating datasets can be accessed via the procedures described [here](https://www.thessgac.org/copy-of-pgi-repository).

## Summary Statistics and PGI Weights
For each phenotype in the Repository, we report GWAS summary statistics and PGI (SBayesR) weights for all SNPs from the largest discovery sample for that analysis, unless the sample includes 23andMe. SNP-level summary statistics from analyses based entirely or in part on 23andMe data can only be reported for up to 10,000 SNPs. Therefore, if the largest GWAS for a phenotype includes 23andMe, we report summary statistics for only the genome-wide significant SNPs from that analysis. In addition, we report summary statistics for all SNPs from the largest GWAS excluding 23andMe. These data will be made available upon publication.

## Measurement-Error-Corrected Estimator
In the initial release of the Repository ([Becker et al. (2021)](https://rdcu.be/cmJnM)), we also propose an approach that improves the interpretability and comparability of research results based on PGIs: to use in place of ordinary least squares (OLS) regression, we derive an estimator that corrects for the  errors-in-variables bias. The estimator produces coefficients in units of the standardized additive SNP factor, which has a more meaningful interpretation than units of some particular PGI. The Python command-line tool implementing the estimator can be found [here](https://github.com/JonJala/pgi_correct).

## Support
This purpose of this repository is to document the code, it is not intended as a tool or library. That being said, we are happy to answer any questions you may have about the code. Before opening an issue, please be sure to read the description of the pipeline in the [paper](https://rdcu.be/cmJnM). 

## Citations
Please include the following citation in any publication based on the Repository PGIs (along with the citations for the GWAS included in the single-trait or multi-trait input GWAS for the PGI) or the measurement error corrected estimator:

Alemu, R. et al. An Updated Polygenic Index Repository: Expanded Phenotypes, New Cohorts, and Improved Causal Inference. bioRxiv 9, 2025.05.14.653986 (2025).
