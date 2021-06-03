#-------------------------------------------------#
# Make additional plots based on the CLEANED QC output
# Rscript --file [FOLDER location of EasyQC results] --ref [path to ref file] 
# Will create supplementary plots for each CLEANED*.gz file found in folder
# Uses reference file with columns: ChrPosID, a1, a2, freq1 
#-------------------------------------------------#

#---------------------------------------------#
# Load necessary packages
#---------------------------------------------#
library("optparse")
library("ggplot2")
library("grid")
library("gridExtra")
library("data.table")

#---------------------------------------------#
# Process user-inputted options 
#---------------------------------------------#
option_list = list(
  make_option(c("--file"),      type="character", default=NULL, help="Folder location of EasyQC results", metavar="character"),
  make_option(c("--ref"),      type="character", default=NULL, help="Path to reference file", metavar="character"),
  make_option(c("--samples"),   type="integer", default=200000, help="Number of samples with which to make plots [default= %default]", metavar="integer"),
  make_option(c("--sdy"),   type="double", default=1.0, help="Standard deviation value of residualized phenotype [default= %default]", metavar="double"),
  make_option(c("--stdx"),   type="logical", default=FALSE, help="TRUE if genotypes are standardized  [default= %default]", metavar="logical"),
  make_option(c("--dom"),   type="logical", default=FALSE, help="TRUE if 23andMe dominance results  [default= %default]", metavar="logical")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Location of results and ref file path must be provided
if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("Must provide path to folder containing EasyQC results!", call.=FALSE)
} else if (is.null(opt$ref)) {
  print_help(opt_parser)
  stop("Must provide path to reference file!", call.=FALSE)
} else {
  cat("Set directory to", opt$file, "\n")
  setwd(opt$file)
}

# Check folder for files matching "CLEANED..gz" to process
d_res_files <- grep("CLEANED.*\\.gz", list.files(), value=TRUE)
if (length(d_res_files) > 0) {
  cat("Found", length(d_res_files), "files to process \n")
  print(d_res_files)
} else {
  stop("No files found that match 'CLEANED...gz'. \n")
}

#---------------------------------------------#
# START
#---------------------------------------------#
cat("-----------------------------------------\n")
cat("QC Step 2 for", opt$file)
Sys.Date()
cat("-----------------------------------------\n")
cat("\n")

# Load Reference Sample 
REF_file <- opt$ref
cat("Loading reference file:", REF_file, "\n")
REF <- fread(REF_file)

#------------------------------------------------------------------------#
# Loop through CLEANED.gz files in the given folder, making QC plots for each:
#------------------------------------------------------------------------#
i <- 1

# Get filename
d_res_file <- d_res_files[i]
in_name  <- tools::file_path_sans_ext(basename(d_res_file), compression = TRUE)
out_name  <- sub("_toQC", "", in_name)

# Log output for each file 
sink(paste0(in_name,".log"), split=TRUE)

# Load data
cat("Loading GWAS data", d_res_file, "\n")
gwas_data <- fread(paste0("gzip -cd ", d_res_file))

# Make plots, sampling 200K SNPs
sample_size <- opt$samples
datasample <- gwas_data[sample(nrow(gwas_data), size=sample_size)]

# TODO: allow a list of sdy's for multiple files with diff var(y)'s? 
sigY <- opt$sdy

# Option for analyzing dominance GWAS results
if (opt$dom == T) {
  datasample[, SE_test := sigY*sqrt(1-(2*EAF*(1-EAF)))/sqrt(2*N*EAF*(1-EAF))]
} else {
    if (opt$stdx == F) {
      datasample[, SE_test := sigY/sqrt(N*2*EAF*(1-EAF))]
    } else {
      # Option for analyzing results for standardized genotype 
      datasample[, SE_test := sigY/sqrt(N)]
    }
}

datasample[, MAF     := ifelse(EAF < 0.5, EAF, 1-EAF)]

#----------------------------------#
# SEN plot (assuming SDY = [user def] and SDX = sqrt(2*maf*(1-maf))
#----------------------------------#
cat("Plotting SEN ... \n")
grid.newpage()

# Run regression to get an estimate of SDy:
datasample[, SE_reg := SE_test/sigY]
fit <- lm(SE ~ 0 + SE_reg, data=datasample)
summary(fit)
sigY_fitted <- coefficients(fit)[[1]]

# Other useful functions 

footnote <- paste0("Number of SNPs sampled: ", sample_size, "\n",
                   "Slope of line assumes SD(y) = ", sigY, "and SD(x) = sqrt(2*maf*(1-maf)). \n",
                   "Estimate of SDy (from regression) = ", sigY_fitted, "\n")

max_x <- max(na.omit(datasample$SE))
max_y <- max(na.omit(datasample$SE_test))

p <- ggplot(data=datasample, aes(x=SE, y=SE_test)) + 
     geom_point(alpha = 0.2) +
     geom_abline(intercept = 0, slope = 1) + 
     theme_bw() + xlab("Actual SE") + ylab("Predicted SE") +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
     xlim(0, max_x) + ylim(0, max_y) + 
     ggtitle(paste0(out_name," Predicted SE vs. Actual SE"))

g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, 
                 gp = gpar(fontface = "italic", fontsize = 12)))

ggsave(plot=g, paste0(out_name,"_SEN.png"), width=6.5, height=6.5)


#----------------------------------#
# Manhattan SE 
#----------------------------------#
cat("Plotting Manhattan SE ... \n")
# Generate actualSE / predicted SE ratio
datasample$SE_ratio <- with(datasample, SE/SE_test)

# extract chr
datasample$CHR <- as.numeric(do.call('rbind',strsplit(datasample$cptid, ":"))[,1])
datasample$POS <- as.numeric(do.call('rbind',strsplit(datasample$cptid, ":"))[,2])

chr1to11 <- datasample[datasample$CHR <= 11 ]
grid.newpage()
footnote <- paste0("Number of SNPs sampled: ", sample_size/2, ". Assume SE = ", sigY, "/sqrt(2*N*maf*(1-maf))")
p <- ggplot(data=chr1to11, aes(x=POS, y=SE_ratio)) + geom_point(alpha=0.4) +
     facet_grid(. ~ CHR) + facet_wrap(~ CHR, ncol = 4, scales="free_x") + 
     theme_bw() + 
     ggtitle(out_name) + xlab("bp") + ylab("SE ratio (Actual/Predicted)")
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 12)))
ggsave(paste0(out_name,"_SEratio_chr1to11.png"), width=6.5, height=6.5, plot=g)

chr12to22 <- datasample[datasample$CHR >= 12 ]
grid.newpage()
footnote <- paste0("Number of SNPs sampled: ", sample_size/2, ". Assume SE = ", sigY, "/sqrt(2*N*maf*(1-maf))")
p <- ggplot(data=chr12to22, aes(x=POS, y=SE_ratio)) + geom_point(alpha=0.4) +
     facet_grid(. ~ CHR) + facet_wrap(~ CHR, ncol = 4, scales="free_x") + 
     theme_bw() + 
     ggtitle(out_name) + xlab("bp") + ylab("SE ratio (Actual/Predicted)")
g <- arrangeGrob(p, bottom = textGrob(footnote, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 12)))
ggsave(paste0(out_name,"_SEratio_chr12to22.png"), width=6.5, height=6.5, plot=g)

sink()
