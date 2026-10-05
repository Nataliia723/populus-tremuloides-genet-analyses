# populus-tremuloides-genet-analyses
R code for genet-level population genetic analyses of trembling aspen 
# Genet-level population genetic analyses of trembling aspen (*Populus tremuloides*)

This repository contains the R code used to perform the genet-level
population genetic analyses reported in the manuscript.

## Analyses

The R script performs the following analyses:

- Construction of consensus genotypes
- Observed heterozygosity (Ho)
- Expected heterozygosity (Hs)
- Inbreeding coefficient (FIS)
- Average number of alleles per locus
- Allelic richness
- Weir & Cockerham FST
- Hedrick's standardized G''ST
- Jost's D

## R script

The main analysis script is:

`genet_level_population_genetics.R`

## Input data

The analysis requires the following input files:

- `second_filters_m4_p80_x0_s3.canonical.unlinked_0.5_100k.vcf`
- `clonal_groups.xlsx`
- `dataset_1.xlsx`

The input datasets are not included in this repository.

## Reproducibility

The R code provided in this repository corresponds to the analyses
reported in the manuscript and supplementary material.

## Code availability

The R scripts used for the analyses are publicly available in this repository.
