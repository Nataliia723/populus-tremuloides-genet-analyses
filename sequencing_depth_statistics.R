# ============================================================
# COVERAGE STATISTICS — CANONICAL SNP DATASET
# Populus tremuloides
#
# Purpose:
# Calculate sequencing-depth statistics reported in the
# Methods in response to Reviewer 3 (R3-17).
#
# Dataset:
# second_filters_m4_p80_x0_s3.canonical.vcf
#
# Expected dataset:
# 86,841 SNPs
# 421 samples
#
# Coverage is calculated from the DP field across all
# non-missing genotype calls.
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGE
# ------------------------------------------------------------

library(vcfR)


# ------------------------------------------------------------
# 2. READ VCF
# ------------------------------------------------------------

vcf <- read.vcfR(
  "C:/Users/Nataliia/Desktop/MAÎTRISE NK/ANALYSES_PEUPLIER/ANALYSES_ROOS/pca-20250624T012301Z-1-001/pca/second_filters_m4_p80_x0_s3.canonical.vcf"
)


# ------------------------------------------------------------
# 3. CHECK DATASET
# ------------------------------------------------------------

n_snps <- nrow(vcf@fix)
n_samples <- ncol(vcf@gt) - 1

cat(
  "Number of SNPs:", n_snps, "\n",
  "Number of samples:", n_samples, "\n"
)

# Stop if a different dataset is accidentally loaded
stopifnot(n_snps == 86841)
stopifnot(n_samples == 421)


# ------------------------------------------------------------
# 4. EXTRACT SEQUENCING DEPTH (DP)
# ------------------------------------------------------------

dp_data <- extract.gt(
  vcf,
  element = "DP",
  as.numeric = TRUE
)

# Convert the DP matrix to a vector
dp_values <- as.vector(dp_data)

# Retain non-missing genotype-depth values
dp_nonmissing <- dp_values[!is.na(dp_values)]


# ------------------------------------------------------------
# 5. COVERAGE STATISTICS
# ------------------------------------------------------------

median_DP <- median(dp_nonmissing)

Q1_DP <- quantile(
  dp_nonmissing,
  probs = 0.25,
  names = FALSE
)

Q3_DP <- quantile(
  dp_nonmissing,
  probs = 0.75,
  names = FALSE
)

mean_DP <- mean(dp_nonmissing)

min_DP <- min(dp_nonmissing)
max_DP <- max(dp_nonmissing)


# ------------------------------------------------------------
# 6. MISSING DP VALUES
# ------------------------------------------------------------

n_DP_total <- length(dp_values)

n_DP_nonmissing <- length(dp_nonmissing)

n_DP_missing <- sum(is.na(dp_values))

percent_DP_missing <- 100 * n_DP_missing / n_DP_total


# ------------------------------------------------------------
# 7. DISPLAY RESULTS
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "CANONICAL SNP DATASET — COVERAGE\n",
  "========================================\n",
  "SNPs:", n_snps, "\n",
  "Samples:", n_samples, "\n",
  "Non-missing DP values:", n_DP_nonmissing, "\n",
  "Missing DP values:", n_DP_missing, "\n",
  "Missing DP (%):", round(percent_DP_missing, 2), "\n\n",
  
  "Median DP:", median_DP, "x\n",
  "Q1:", Q1_DP, "x\n",
  "Q3:", Q3_DP, "x\n",
  "IQR:", Q1_DP, "-", Q3_DP, "x\n",
  "Mean DP:", round(mean_DP, 2), "x\n",
  "Minimum DP:", min_DP, "x\n",
  "Maximum DP:", max_DP, "x\n"
)


# ------------------------------------------------------------
# 8. ADDITIONAL QUANTILES
# ------------------------------------------------------------

DP_quantiles <- quantile(
  dp_nonmissing,
  probs = c(
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99
  )
)

print(DP_quantiles)


# ------------------------------------------------------------
# 9. RESULTS REPORTED IN THE MANUSCRIPT
# ------------------------------------------------------------

cat(
  "\nMANUSCRIPT VALUE:\n",
  "Median sequencing depth =",
  median_DP,
  "x (IQR:",
  Q1_DP,
  "-",
  Q3_DP,
  "x)\n"
)


# ------------------------------------------------------------
# 10. SESSION INFORMATION
# ------------------------------------------------------------

sessionInfo()