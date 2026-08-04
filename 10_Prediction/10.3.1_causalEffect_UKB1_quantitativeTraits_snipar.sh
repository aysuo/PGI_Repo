#!/bin/bash
# ------------------------------------------------------------------------------
# Script: 10.3.run_causalEffect_UKB1_quantitativeTraits_snipar.sh
# ------------------------------------------------------------------------------
# Description:
# This script estimates causal (direct) effects of polygenic indexes (PGIs) for 
# quantitative traits in the UKB1 cohort using the SNIPAR package 
# (Young et al., 2022). It runs SNIPAR separately for each phenotype listed 
# and compiles proband effect estimates into one file.
#
# Associated manuscript:
# Alemu et al. (2025) "An Updated Polygenic Index Repository: Expanded Phenotypes,
# New Cohorts, and Improved Causal Inference."
# ------------------------------------------------------------------------------

# --- 1. Source paths and environment variables ---
source paths_PGIrepo_withinFam   # This defines PGI_RepoV2 and related directories
source $snipar_venv              # Activate SNIPAR virtual environment

# --- 2. Define core paths ---
PYTHON_BIN="$PGI_RepoV2/doc/python/bin/python3.9"      # Python interpreter in venv
phenotype_list_file="$PGI_RepoV2/scripts/7_withinFamily_related/ukb1_pheno_list"
output_directory="$PGI_RepoV2/processed/sumStats/within_family_related/UKB1"

covariate_file="$UKB1_covariates"
phenotype_dir="$UKB1_phenos"
pgi_dir="$UKB1_pgiparents"

# --- 3. Loop through phenotypes and run SNIPAR ---
while IFS= read -r pheno
do
    echo "Processing phenotype: $pheno"

    # Define phenotype and PGI file paths
    phenotype_file="${phenotype_dir}/${pheno}.pheno"
    pgi_file="${pgi_dir}/PGS_UKB1_parental_${pheno}-single_SBayesR.pgs.txt"

    # Temporary phenotype file (replace missing with NA)
    temp_phenotype_file="temp_${pheno}.pheno"
    awk '{if (NF == 2) print $1, $2, "NA"; else print $0}' "$phenotype_file" > "$temp_phenotype_file"

    # Temporary covariate file (keep FID, IID, and PCs)
    temp_covariate_file="temp_covar_${pheno}.txt"
    awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22}' \
        "$covariate_file" > "$temp_covariate_file"

    # Run SNIPAR (Young et al., 2022) for causal effect estimates
    $PYTHON_BIN $pgs_py direct \
        --pgs "$pgi_file" \
        --phenofile "$temp_phenotype_file" \
        --covar "$temp_covariate_file" \
        --scale_pgs --scale_phen

    # Move outputs to the appropriate directory
    for suffix in "1.effects.txt" "2.effects.txt" "1.vcov.txt" "2.vcov.txt"
    do
        if [ -f direct.$suffix ]; then
            mv direct.$suffix ${output_directory}/${pheno}_direct.$suffix
        fi
    done

    # Clean up temporary files
    rm "$temp_phenotype_file" "$temp_covariate_file"

    echo "Completed processing for $pheno. Outputs saved in ${output_directory}"
done < "$phenotype_list_file"

echo "All phenotypes processed."

# --- 4. Compile proband effects into one file (Python block) ---
$PYTHON_BIN << EOF
import os
import pandas as pd

# Define paths
pgi_repo = os.environ.get("PGI_RepoV2")
effect_files_dir = os.path.join(pgi_repo, "processed/sumStats/within_family_related/UKB1")
phenotype_list_file = os.path.join(pgi_repo, "scripts/7_withinFamily_related/ukb1_pheno_list")
output_file = os.path.join(effect_files_dir, "compiled_proband_effects.txt")

# Read phenotypes
with open(phenotype_list_file, 'r') as f:
    phenotypes = f.read().splitlines()

compiled_data = []
for pheno in phenotypes:
    effect_file = os.path.join(effect_files_dir, f"{pheno}_direct.2.effects.txt")
    if os.path.exists(effect_file):
        with open(effect_file, 'r') as ef:
            for line in ef:
                if line.startswith("proband"):
                    parts = line.split()
                    compiled_data.append([pheno, parts[1], parts[2]]) # phenotype, coefficient, SE
                    break

# Save compiled table
df = pd.DataFrame(compiled_data, columns=["Phenotype","Coefficient","StandardError"])
df.to_csv(output_file, sep="\t", index=False)
print(f"Compiled proband effects saved to {output_file}")
EOF

