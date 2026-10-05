# ============================================================
# GENET-LEVEL POPULATION GENETIC ANALYSES
# Trembling aspen (Populus tremuloides)
# ============================================================
#
# Analyses performed at the genet level to avoid
# pseudoreplication caused by repeated sampling of ramets
# belonging to the same genet.
#
# Dataset:
#   - 421 genotyped samples
#   - 25 multi-ramet genets (413 samples)
#   - 8 samples not assigned to a multi-ramet genet
#   - 33 distinct genotypes used for population analyses
#   - 48,976 unlinked SNPs
#
# Genet-level sample sizes:
#   Berry    = 6
#   Desboues = 2
#   Rapide   = 25
#
# Analyses:
#   1. Construction of consensus genotypes
#   2. Genetic diversity: Ho, Hs, FIS
#   3. Average number of alleles per locus
#   4. Allelic richness
#   5. Weir & Cockerham FST
#   6. Hedrick's standardized G''ST
#   7. Jost's D
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

library(vcfR)
library(openxlsx)
library(dplyr)
library(adegenet)
library(hierfstat)
library(mmod)


# ------------------------------------------------------------
# 2. INPUT FILES
# ------------------------------------------------------------

# Unlinked SNP dataset used for population genetic analyses
vcf_unlinked <- read.vcfR(
  "second_filters_m4_p80_x0_s3.canonical.unlinked_0.5_100k.vcf"
)

# Assignment of samples to the 25 multi-ramet genets
clone_groups <- read.xlsx("clonal_groups.xlsx")

# Sample metadata
dataset_1 <- read.xlsx("dataset_1.xlsx")


# ------------------------------------------------------------
# 3. EXTRACT GENOTYPES
# ------------------------------------------------------------

gt_unlinked <- extract.gt(
  vcf_unlinked,
  element = "GT",
  as.numeric = FALSE
)

dim(gt_unlinked)
# Expected: 48,976 SNPs x 421 samples


# ------------------------------------------------------------
# 4. IDENTIFY SAMPLES NOT ASSIGNED TO MULTI-RAMET GENETS
# ------------------------------------------------------------

vcf_samples <- colnames(gt_unlinked)

unique_samples <- setdiff(
  vcf_samples,
  clone_groups$Sample
)

unique_samples

# Expected:
# ptrem_014
# ptrem_080
# ptrem_174
# ptrem_293
# ptrem_349
# ptrem_410
# ptrem_539
# ptrem_540

length(unique_samples)
# Expected: 8


# ------------------------------------------------------------
# 5. CONSTRUCT CONSENSUS GENOTYPES FOR THE 25 MULTI-RAMET GENETS
# ------------------------------------------------------------

# Function returning the most frequent non-missing genotype.
# In case of a tie, the first modal genotype is retained.

consensus_fun <- function(x) {
  
  x <- x[
    !is.na(x) &
      x != "./." &
      x != ".|." &
      x != "."
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  tab <- table(x)
  
  names(tab)[which.max(tab)]
}


genet_ids <- sort(unique(clone_groups$Clone_Group))

consensus_unlinked <- matrix(
  NA_character_,
  nrow = nrow(gt_unlinked),
  ncol = length(genet_ids)
)

rownames(consensus_unlinked) <- rownames(gt_unlinked)
colnames(consensus_unlinked) <- paste0("G", genet_ids)


for (i in seq_along(genet_ids)) {
  
  g <- genet_ids[i]
  
  samples_g <- clone_groups$Sample[
    clone_groups$Clone_Group == g
  ]
  
  samples_g <- intersect(
    samples_g,
    colnames(gt_unlinked)
  )
  
  consensus_unlinked[, i] <- apply(
    gt_unlinked[, samples_g, drop = FALSE],
    1,
    consensus_fun
  )
}


# ------------------------------------------------------------
# 6. ADD THE EIGHT UNIQUE GENOTYPES
# ------------------------------------------------------------

unique_unlinked <- gt_unlinked[
  ,
  unique_samples,
  drop = FALSE
]

colnames(unique_unlinked) <- paste0(
  "U_",
  sub("ptrem_", "", unique_samples)
)


gt_33_unlinked <- cbind(
  consensus_unlinked,
  unique_unlinked
)

dim(gt_33_unlinked)
# Expected: 48,976 x 33


# ------------------------------------------------------------
# 7. CREATE GENET-LEVEL METADATA
# ------------------------------------------------------------

# Site of each multi-ramet genet
genet_site <- clone_groups %>%
  left_join(
    dataset_1[, c("Sample", "Site")],
    by = "Sample"
  ) %>%
  group_by(Clone_Group) %>%
  summarise(
    Site = first(Site),
    N_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Genotype = paste0("G", Clone_Group),
    Type = "Multi-ramet genet"
  )


# Metadata for the eight unique genotypes
unique_site <- dataset_1 %>%
  filter(Sample %in% unique_samples) %>%
  transmute(
    Genotype = paste0(
      "U_",
      sub("ptrem_", "", Sample)
    ),
    Site = Site,
    Type = "Unique genotype",
    N_samples = 1
  )


metadata_33 <- bind_rows(
  genet_site %>%
    select(Genotype, Site, Type, N_samples),
  unique_site
)


# Reorder metadata to exactly match genotype matrix
metadata_33_unlinked <- metadata_33[
  match(
    colnames(gt_33_unlinked),
    metadata_33$Genotype
  ),
]


# Verification
stopifnot(
  all(
    metadata_33_unlinked$Genotype ==
      colnames(gt_33_unlinked)
  )
)

table(metadata_33_unlinked$Site)

# Expected:
# Berry     6
# Desboues  2
# Rapide   25


# ------------------------------------------------------------
# 8. CONVERT GENOTYPES FOR HIERFSTAT
# ------------------------------------------------------------

# hierfstat diploid genotype coding:
# 0/0 -> 11
# 0/1 -> 12
# 1/0 -> 12
# 1/1 -> 22

convert_hf <- function(x) {
  
  x[x %in% c("./.", ".|.", ".")] <- NA
  
  x[x %in% c("0/0", "0|0")] <- "11"
  
  x[x %in% c(
    "0/1", "1/0",
    "0|1", "1|0"
  )] <- "12"
  
  x[x %in% c("1/1", "1|1")] <- "22"
  
  as.numeric(x)
}


gt_hf_numeric <- apply(
  gt_33_unlinked,
  2,
  convert_hf
)

# Transpose:
# individuals/genets in rows
# SNPs in columns

gt_hf_numeric <- t(gt_hf_numeric)


# ------------------------------------------------------------
# 9. ASSIGN POPULATIONS
# ------------------------------------------------------------

site_factor <- factor(
  metadata_33_unlinked$Site,
  levels = c(
    "Berry",
    "Desboues",
    "Rapide"
  )
)

pop_code <- as.integer(site_factor)


hf_33 <- data.frame(
  pop = pop_code,
  gt_hf_numeric,
  check.names = FALSE
)

rownames(hf_33) <-
  metadata_33_unlinked$Genotype


dim(hf_33)
# Expected: 33 x 48,977


# ============================================================
# 10. GENETIC DIVERSITY
# ============================================================

basic_33 <- basic.stats(hf_33)


# Mean observed heterozygosity
Ho_33 <- colMeans(
  basic_33$Ho,
  na.rm = TRUE
)

# Mean expected heterozygosity
Hs_33 <- colMeans(
  basic_33$Hs,
  na.rm = TRUE
)

# Mean FIS across loci
Fis_33 <- colMeans(
  basic_33$Fis,
  na.rm = TRUE
)


diversity_33 <- data.frame(
  Site = c(
    "Berry",
    "Desboues",
    "Rapide"
  ),
  N_genets = c(6, 2, 25),
  Ho = Ho_33,
  Hs = Hs_33,
  FIS = Fis_33
)

diversity_33


# Expected approximately:
#
# Site       N_genets    Ho       Hs       FIS
# Berry          6     0.1526   0.1733    0.0734
# Desboues       2     0.1176   0.0564   -0.4964
# Rapide        25     0.1523   0.1720    0.0752


# ============================================================
# 11. AVERAGE NUMBER OF ALLELES PER LOCUS
# ============================================================

# basic.stats() stores allele frequencies as a list,
# with one element per SNP.

n_alleles_locus <- t(
  sapply(
    basic_33$pop.freq,
    function(x) {
      
      apply(
        x,
        2,
        function(freq) {
          sum(freq > 0, na.rm = TRUE)
        }
      )
    }
  )
)


# Check that no population-locus combination
# was counted as having zero observed alleles

zero_allele_counts <- colSums(
  n_alleles_locus == 0,
  na.rm = TRUE
)

zero_allele_counts
# Expected: 0 0 0


mean_alleles_33 <- colMeans(
  n_alleles_locus,
  na.rm = TRUE
)

mean_alleles_33

# Expected:
# Berry     1.559886
# Desboues  1.139762
# Rapide    1.925780


diversity_33$Mean_alleles_per_locus <-
  mean_alleles_33


# Final Table S6 values
diversity_33


# ============================================================
# 12. ALLELIC RICHNESS
# ============================================================

AR_33 <- allelic.richness(hf_33)

Ar_mean_33 <- colMeans(
  AR_33$Ar,
  na.rm = TRUE
)

Ar_mean_33

# Expected:
# Berry     1.171107
# Desboues  1.103626
# Rapide    1.171515


AR_33$min.all
# Expected: 2

# NOTE:
# Allelic richness is rarefied to the minimum sample size,
# which is only two genets because Desboues contains n = 2.
# Therefore, these estimates should be interpreted cautiously.
# Allelic richness was not substituted for the "average number
# of alleles per locus" reported in Table S6.


# ============================================================
# 13. WEIR & COCKERHAM FST
# ============================================================

fst_wc_33 <- wc(hf_33)

fst_wc_33$FST
# Expected: 0.02144234

fst_wc_33$FIS
# Expected: 0.1128775


# Pairwise Weir & Cockerham FST

fst_pairwise_33 <- pairwise.WCfst(hf_33)

round(
  fst_pairwise_33,
  4
)

# Expected:
#
#             Berry  Desboues  Rapide
# Berry          NA   0.0839   0.0026
# Desboues   0.0839       NA   0.0639
# Rapide     0.0026   0.0639       NA


# ============================================================
# 14. CONVERT TO GENIND FOR COMPLEMENTARY
#     DIFFERENTIATION METRICS
# ============================================================

# Build a temporary VCF containing only the 33
# genet-level genotype columns.

vcf_33 <- vcf_unlinked

vcf_33@gt <- cbind(
  FORMAT = vcf_unlinked@gt[, "FORMAT"],
  gt_33_unlinked
)

genind_33 <- vcfR2genind(vcf_33)

indNames(genind_33) <-
  metadata_33_unlinked$Genotype

pop(genind_33) <-
  factor(metadata_33_unlinked$Site)


# Verification
nInd(genind_33)
# Expected: 33

nLoc(genind_33)
# Expected: 48,976

table(pop(genind_33))
# Expected:
# Berry     6
# Desboues  2
# Rapide   25


# ============================================================
# 15. HEDRICK'S STANDARDIZED G''ST
# ============================================================

Gst_hedrick_33 <-
  Gst_Hedrick(genind_33)

Gst_hedrick_33

# Expected global G''ST:
# 0.09103191


Gst_hedrick_pairwise_33 <-
  pairwise_Gst_Hedrick(genind_33)

round(
  Gst_hedrick_pairwise_33,
  4
)

# Expected:
#
#                   Berry  Desboues
# Desboues         0.1482
# Rapide           0.0114   0.1312


# ============================================================
# 16. JOST'S D
# ============================================================

D_jost_33 <- D_Jost(genind_33)

# Global Jost's D
D_jost_33$global.het

# Expected:
# 0.0149 approximately


# Pairwise Jost's D
D_jost_pairwise_33 <-
  pairwise_D(genind_33)

round(
  D_jost_pairwise_33,
  4
)

# Expected:
#
#                   Berry  Desboues
# Desboues         0.0237
# Rapide           0.0020   0.0209


# ============================================================
# 17. FINAL DIFFERENTIATION SUMMARY
# ============================================================

differentiation_summary <- data.frame(
  
  Metric = c(
    "Weir & Cockerham FST",
    "Hedrick standardized G''ST",
    "Jost D"
  ),
  
  Overall = c(
    fst_wc_33$FST,
    Gst_hedrick_33,
    D_jost_33$global.het
  ),
  
  Berry_Desboues = c(
    fst_pairwise_33[1, 2],
    Gst_hedrick_pairwise_33["Desboues", "Berry"],
    D_jost_pairwise_33["Desboues", "Berry"]
  ),
  
  Berry_Rapide = c(
    fst_pairwise_33[1, 3],
    Gst_hedrick_pairwise_33["Rapide", "Berry"],
    D_jost_pairwise_33["Rapide", "Berry"]
  ),
  
  Desboues_Rapide = c(
    fst_pairwise_33[2, 3],
    Gst_hedrick_pairwise_33["Rapide", "Desboues"],
    D_jost_pairwise_33["Rapide", "Desboues"]
  )
)

differentiation_summary


# ============================================================
# 18. EXPORT TABLES
# ============================================================

write.csv(
  diversity_33,
  "Table_S6_genet_level_genetic_diversity.csv",
  row.names = FALSE
)

write.csv(
  differentiation_summary,
  "Genet_level_genetic_differentiation.csv",
  row.names = FALSE
)

write.csv(
  metadata_33_unlinked,
  "Genet_level_sample_metadata_33_genotypes.csv",
  row.names = FALSE
)


# ============================================================
# END
# ============================================================
