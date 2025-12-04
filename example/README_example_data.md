# PALM-mbQTL Example Dataset (SeqDepth = sum of ABD)

This toy dataset is designed to let users run the entire PALM-mbQTL
pipeline. SeqDepth is computed **as the sum of all abundance features**
for each sample, which matches typical microbiome workflows.

## Files
- abd.tsv — microbiome abundances
- cov.tsv — covariates + SeqDepth = g_Blautia + g_Roseburia
- pheno_all_taxa.txt — merged phenotype file
- example.ped / example.map — tiny PLINK dataset

## Generate BED/BIM/FAM
```
plink --ped example.ped --map example.map --make-bed --out example
```

## Convert to GDS (optional)
```
library(SNPRelate)
snpgdsPED2GDS("example.ped", "example.map", "example.gds")
```

