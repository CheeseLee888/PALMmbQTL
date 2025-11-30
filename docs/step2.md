---
title: Step 2 – Score Tests
nav_order: 7
---

<!--
Detailed description of Step 2: genome-wide score test computation.
-->

# Step 2: Genome-wide Score Tests

Step 2 computes genome-wide score statistics and p-values for a single
microbiome feature using the null model fitted in Step 1.

## Script location

The main script for Step 2 is:

```text
extdata/step2_scoreTest.R
```

View the available options:

```bash
cd PALMmbQTL/extdata
pixi run --manifest-path=../pixi.toml Rscript step2_scoreTest.R --help
```

## Typical arguments

Exact names may differ; check `--help` for the current version.
Common options:

- `--inFile`          : genotype data (GDS file)
- `--GMMATmodelFile`  : RDS file from Step 1 (`*_step1.rda`)
- `--SAIGEOutputFile` : output text file with SNP-level results
- `--LOCO`            : logical flag for leave-one-chromosome-out (LOCO) analysis

Additional options may control variant filtering, chromosome selection,
or parallelization.

## Example (pixi)

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step2_scoreTest.R   --inFile=/path/to/mydata/stool_bialleic_merged_data.gds   --GMMATmodelFile=/path/to/mydata/output/g_Blautia_step1.rda   --SAIGEOutputFile=/path/to/mydata/output/g_Blautia_allchr_step2.txt   --LOCO=FALSE
```

## Example (Docker)

```bash
docker run --rm -v /path/to/mydata:/data palmmbqtl:latest   step2_scoreTest.R   --inFile=/data/stool_bialleic_merged_data.gds   --GMMATmodelFile=/data/output/g_Blautia_step1.rda   --SAIGEOutputFile=/data/output/g_Blautia_allchr_step2.txt   --LOCO=FALSE
```

## Output

The output file (`*_allchr_step2.txt`) is typically a tab-delimited text file
with one row per variant. Columns may include:

- Chromosome and position
- Variant ID
- Allele information
- Score statistic and variance
- p-value

Use your preferred downstream tools to visualize or filter these results.

This concludes the core PALM-mbQTL pipeline. You can now:

- Repeat Steps 1–2 for additional microbiome features.
- Combine results across features or cohorts as needed.
