---
title: Step 1 – Null Model
nav_order: 6
---

<!--
Detailed description of Step 1: fitting the null Poisson mixed model.
-->

# Step 1: Fitting the Null Poisson Mixed Model

Step 1 fits a **Poisson mixed null model** for a single microbiome feature,
such as `g_Blautia`. The fitted model is saved and reused for genome-wide
score tests in Step 2.

## Script location

The main script for Step 1 is:

```text
extdata/step1_fitNULL.R
```

To view the command-line options:

```bash
cd PALMmbQTL/extdata
pixi run --manifest-path=../pixi.toml Rscript step1_fitNULL.R --help
```

## Typical arguments

Common command-line options include:

- `--phenoFile` : phenotype file (tab-delimited)
- `--sampleCol` : name of the sample ID column
- `--yCol`      : name of the microbiome feature column (e.g. `g_Blautia`)
- `--covarCols` : comma-separated list of covariate column names
- `--grmFile`   : GRM RDS file (optional; identity is used if omitted)
- `--outFile`   : RDS file to save the fitted null model

## Example (pixi)

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step1_fitNULL.R   --phenoFile=/path/to/mydata/pheno_all_taxa.txt   --sampleCol=sample_id   --yCol=g_Blautia   --covarCols=age,sex,PC1,PC2   --grmFile=/path/to/mydata/output/stool_GRM.rds   --outFile=/path/to/mydata/output/g_Blautia_step1.rda
```

## Example (Docker)

```bash
docker run --rm -v /path/to/mydata:/data palmmbqtl:latest   step1_fitNULL.R   --phenoFile=/data/pheno_all_taxa.txt   --sampleCol=sample_id   --yCol=g_Blautia   --covarCols=age,sex,PC1,PC2   --grmFile=/data/output/stool_GRM.rds   --outFile=/data/output/g_Blautia_step1.rda
```

## Output

The output file (`*_step1.rda`) contains the fitted null model object.
Typically this includes:

- Fixed effect coefficients
- Variance component estimates
- Information needed to compute score statistics in Step 2

Next: [Step 2 – Score Tests](step2.md).
