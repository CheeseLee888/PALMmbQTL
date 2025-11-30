---
title: Overview
nav_order: 2
---

<!--
High-level conceptual overview of the PALM-mbQTL analysis pipeline.
Adapt this page when your methodology or implementation changes.
-->

# PALM-mbQTL Analysis Overview

PALM-mbQTL provides a multi-step workflow for microbiome QTL mapping:

1. **Prepare input data**
   - Genotype (GDS or PLINK)
   - Phenotype (microbiome feature table + covariates)
   - Optional GRM

2. **Step 0 (optional but recommended)**  
   Construct a genetic relationship matrix (GRM) from genotype data.

3. **Step 1**  
   Fit a Poisson mixed null model for one microbiome feature.

4. **Step 2**  
   Use the fitted null model to compute genome-wide score statistics.

## Model components

### Genotype and GRM

- Genotype data are stored either as:
  - SeqArray/SNP GDS (recommended), or
  - PLINK bed/bim/fam, which can be converted to GDS.
- The GRM summarizes sample relatedness and is used as a random effect.

### Null Poisson mixed model (Step 1)

For a single microbiome feature (e.g. genus `g_Blautia`), PALM-mbQTL fits:

- Fixed effects: covariates such as age, sex, PCs, batch, etc.
- Random effects:
  - Feature-shared random intercept to account for compositional structure.
  - Optional GRM-based random effect to account for sample relatedness.

The fitted model is saved for reuse in Step 2.

### Score tests (Step 2)

Using the fitted null model, PALM-mbQTL computes score statistics and
p-values for each variant. This approach avoids re-fitting the full model
for every SNP and scales to large genotype matrices.

## Input files

You typically need:

- **Phenotype file**

  - One row per sample.
  - Columns for sample ID, microbiome features, and covariates.

- **Genotype file**

  - GDS file, or
  - PLINK bed/bim/fam (converted internally to GDS if supported).

- **GRM file** (optional)

  - RDS file containing a genetic relationship matrix and sample IDs.
  - The set and order of sample IDs must match the phenotype file.

## Computational considerations

- Step 0 (GRM) is run once per cohort / genotype matrix.
- Step 1 is run once per microbiome feature and can be parallelized.
- Step 2 is fast and can be parallelized by chromosome or variant blocks.

Next steps:

- [Installation](installation.md)
- [Quick Start guide](quickstart.md)
