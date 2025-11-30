---
title: Quick Start
nav_order: 4
---

<!--
This page demonstrates a minimal end-to-end run for one microbiome feature.
Adapt file paths and argument names to your actual implementation.
-->

# Quick Start

This Quick Start example shows how to run PALM-mbQTL for a single
microbiome feature, assumed to be named `g_Blautia`.

The pipeline consists of:

1. Step 0: generate a GRM (optional but recommended)
2. Step 1: fit the null Poisson mixed model
3. Step 2: compute genome-wide score statistics

Throughout this page we assume that:

- The PALMmbQTL repository is cloned at `/path/to/PALMmbQTL`.
- Input files are under `/path/to/mydata`.
- Output files are written back to `/path/to/mydata/output`.

## 1. Prepare data

You need:

- A genotype file (GDS or PLINK).
- A phenotype file with:
  - sample ID column,
  - a column for `g_Blautia`, and
  - covariates such as age, sex, PCs.

Example phenotype table:

```text
sample_id   g_Blautia   age  sex  PC1   PC2
S1          12          35   M    0.12 -0.03
S2          20          42   F   -0.02  0.11
S3          18          29   F    0.05  0.00
```

## 2. Step 0 – Generate a GRM

Using pixi:

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step0_generateGRM.R   --genoFile=/path/to/mydata/stool_bialleic_merged_data.gds   --grmFile=/path/to/mydata/output/stool_GRM.rds
```

Using Docker:

```bash
docker run --rm -v /path/to/mydata:/data palmmbqtl:latest   step0_generateGRM.R   --genoFile=/data/stool_bialleic_merged_data.gds   --grmFile=/data/output/stool_GRM.rds
```

If you omit `--grmFile`, or if no GRM is available, the identity matrix
may be used as a default depending on your implementation.

## 3. Step 1 – Fit null model for `g_Blautia`

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step1_fitNULL.R   --phenoFile=/path/to/mydata/pheno_all_taxa.txt   --sampleCol=sample_id   --yCol=g_Blautia   --covarCols=age,sex,PC1,PC2   --grmFile=/path/to/mydata/output/stool_GRM.rds   --outFile=/path/to/mydata/output/g_Blautia_step1.rda
```

The file `g_Blautia_step1.rda` stores the fitted null model and is required
for Step 2.

## 4. Step 2 – Genome-wide score tests

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step2_scoreTest.R   --inFile=/path/to/mydata/stool_bialleic_merged_data.gds   --GMMATmodelFile=/path/to/mydata/output/g_Blautia_step1.rda   --SAIGEOutputFile=/path/to/mydata/output/g_Blautia_allchr_step2.txt   --LOCO=FALSE
```

After completion, the file `g_Blautia_allchr_step2.txt` should contain
per-SNP statistics and p-values.

For more detailed explanations of each step, see:

- [Step 0 – GRM Construction](step0.md)
- [Step 1 – Null Model Fitting](step1.md)
- [Step 2 – Score Tests](step2.md)
