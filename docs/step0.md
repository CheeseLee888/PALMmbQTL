---
title: Step 0 – GRM Construction
nav_order: 5
---

<!--
Detailed documentation for GRM construction.
You can adjust naming/arguments to match your actual script.
-->

# Step 0: Constructing a Genetic Relationship Matrix (GRM)

Step 0 generates a genetic relationship matrix (GRM) from genotype data.
This step is optional but recommended when individuals may be related
or when population structure needs to be modeled.

## Script location

The helper script is typically located at:

```text
extdata/step0_generateGRM.R
```

You can view its usage via:

```bash
cd PALMmbQTL/extdata
pixi run --manifest-path=../pixi.toml Rscript step0_generateGRM.R --help
```

## Typical arguments

Exact argument names may differ; consult `--help` for the version you use.
Common options include:

- `--genoFile` : path to genotype data (GDS or PLINK prefix)
- `--grmFile`  : output RDS file containing the GRM
- `--maf`      : minimum minor allele frequency filter (optional)
- `--missing`  : missingness threshold (optional)

## Example (GDS input)

```bash
cd /path/to/PALMmbQTL/extdata

pixi run --manifest-path=../pixi.toml   Rscript step0_generateGRM.R   --genoFile=/path/to/mydata/stool_bialleic_merged_data.gds   --grmFile=/path/to/mydata/output/stool_GRM.rds
```

## Output format

The GRM RDS file is expected to contain at least:

- A vector of sample IDs (e.g. `sample.id`)
- A square matrix `K` representing pairwise relationships

The order of sample IDs must match the ordering used in the phenotype file
for Step 1. The PALM-mbQTL implementation checks that the IDs are consistent
and will throw an error if they are not.

Once the GRM is constructed, proceed to
[Step 1 – Null Model Fitting](step1.md).
