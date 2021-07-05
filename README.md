# PGI_Repo
This repository contains the complete pipeline for the creation of polygenic indexes (PGIs) for the SSGAC Polygenic Index Repository, an initiative that makes PGIs for a wide range of traits available for a number of datasets that may be useful to social scientists. The Repository currently contains single- and/or multi-trait (MTAG) PGIs for [47 phenotypes](https://www.dropbox.com/s/bo68hgfp8d7ibgk/Phenotypes.pdf?dl=0) in [11 datasets](https://www.dropbox.com/s/psrqvr5qoop4zye/Datasets.pdf?dl=0). To maximize prediction accuracy of the PGIs, we meta-analysed summary statistics from multiple sources, including several novel large-scale GWASs conducted in UK Biobank and the personal genomics company 23andMe. Therefore, almost all PGIs in our initial release perform at least as well as currently available PGIs in terms of prediction accuracy. Please see [Becker et al. (2021)](https://rdcu.be/cmJnM)  and the [User Guide](https://8a649e3c-e96e-4fce-86bb-b117a3d1270f.filesusr.com/ugd/2f9665_2d32cc8c8b63449d994daefbdd6c77b1.pdf) for a detailed description of the pipeline. 

The Repository will be updated regularly with additional PGIs and datasets. If you are interested in participating in the Repository, please reach out to contact@ssgac.org. 

## Frequently Asked Questions (FAQs)
For a less technical description of the paper and of how PGIs should—and should not—be interpreted and used, see these [frequently asked questions](https://www.thessgac.org/faqs).

## PGI Access Procedures
PGIs in the participating datasets can be accessed via the procedures described [here](https://www.thessgac.org/copy-of-pgi-repository).

## Summary Statistics and PGI Weights
For each phenotype in the Repository, we report GWAS and MTAG summary statistics and PGI (LDpred) weights for all SNPs from the largest discovery sample for that analysis, unless the sample includes 23andMe. SNP-level summary statistics from analyses based entirely or in part on 23andMe data can only be reported for up to 10,000 SNPs. Therefore, if the largest GWAS or MTAG analysis for a phenotype includes 23andMe, we report summary statistics for only the genome-wide significant SNPs from that analysis. In addition, we report summary statistics for all SNPs from the largest GWAS or MTAG analysis excluding 23andMe. These data will be made available upon publication.

## Measurement-Error-Corrected Estimator
In [Becker et al. (2021)](https://rdcu.be/cmJnM), we also propose an approach that improves the interpretability and comparability of research results based on PGIs: to use in place of ordinary least squares (OLS) regression, we derive an estimator that corrects for the  errors-in-variables bias. The estimator produces coefficients in units of the standardized additive SNP factor, which has a more meaningful interpretation than units of some particular PGI. The Python command-line tool implementing the estimator can be found [here](https://github.com/JonJala/pgi_correct).

## Support
This purpose of this repository is to document the code, it is not intended as a tool or library. That being said, we are happy to answer any questions you may have about the code. Before opening an issue, please be sure to read the description of the pipeline in the [paper](https://rdcu.be/cmJnM). 

## Citations
Please include the following citation in any publication based on the Repository PGIs (along with the citations for the GWAS included in the single-trait or multi-trait input GWAS for the PGI) or the measurement error corrected estimator:

Becker, J., Burik, C.A.P., Goldman, G., Wang, N., Jayashankar, H., Bennett, M., Belsky, D.W., Karlsson Linnér, R., Ahlskog, R., Kleinman, A., Hinds, D.A., 23andMe Research Group, Caspi, A., Corcoran, D.L., Moffitt, T.E., Poulton, R., Sugden, K., Williams, B.S., Harris, K.M., Steptoe, A., Ajnakina, O., Milani, L., Esko, T., Iacono, W.G., McGue, T., Magnusson, P.K.E., Mallard, T.T., Harden, K.P., Tucker-Drob, E.M., Herd, P., Freese, J., Young, A., Beauchamp, J.P., Koellinger, P.D., Oskarsson, S., Johannesson, M., Visscher, P.M., Meyer, M.N., Laibson, D., Cesarini, D., Benjamin, D.J., Turley, P., and Okbay, A. (2021). Resource Profile and User Guide of the Polygenic Index Repository. *Nature Human Behaviour*. Published online June 17. doi:10.1038/s41562-021-01119-3.
