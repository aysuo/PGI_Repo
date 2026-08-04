
# load libraries
packages <- c("rhdf5")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# Parse arguments
args=commandArgs(trailingOnly=TRUE)
cohort=args[1]
gfDir=args[2]
outDir=args[3]


bimmerged <- data.frame()

if ( cohort %in% c("ALSPAC", "PSID") ) {
    cols <- c("chr","snpid","yyy","pos","A1","A2")
} else if ( cohort %in% c("AH", "UKB1") ) {
    cols <- c("chr","snpid","pos","A1","A2")
} else {
    cols <- c("snpid","alternative_snpid","chr","pos","n_alleles","alleles","xxx")
}

for (i in 1:22) {

    if ( cohort == "UKB1" ) {
        hdf5name <- paste0(gfDir,"/chr_", i, ".hdf5")
    } else {
        hdf5name <- paste0(gfDir,"/",cohort,"_parental_chr", i, ".hdf5")
    }

    bim <- h5read(hdf5name, "/bim_values")
    bim <- as.data.frame(t(bim))
    bimmerged <- rbind(bim, bimmerged)
    colnames(bim) <- cols
    write.table(bim,paste0(outDir,"/chr",i,".bim"), quote=F, row.names=F)
}

colnames(bimmerged) <- cols
write.table(bimmerged,paste0(outDir,"/",cohort,".bim"), quote=F, row.names=F)