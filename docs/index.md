---
title: PALM-mbQTL
nav_order: 1
---

<!--
This is the home page of the PALM-mbQTL documentation.
Edit the text below to customize the landing page for your project.
-->

# PALM-mbQTL

PALM-mbQTL is an R-based framework (package name: `SAIGEQTL`) for
**fast and scalable Poisson mixed–model microbiome QTL (mbQTL) mapping**.

The method:

- Models microbiome features (taxa, genes, pathways, etc.) as discrete counts.
- Uses a feature-shared random intercept to account for compositional correlation.
- Incorporates a genetic relationship matrix (GRM) to handle related samples.
- Reuses fitted null models to compute genome-wide score statistics efficiently.

## Getting started

1. Read the [Overview](overview.md) to understand the analysis workflow.
2. Install the package and environment: see [Installation](installation.md).
3. Follow the [Quick Start](quickstart.md) guide to run a minimal example.
4. Dive into individual steps:
   - [Step 0 – GRM Construction](step0.md)
   - [Step 1 – Null Model Fitting](step1.md)
   - [Step 2 – Genome-wide Score Tests](step2.md)

## Citation

> Li P, Zhou W, Tang Z, *et al.*  
> **PALM-mbQTL: Poisson mixed–model framework for microbiome QTL mapping.**  
> (Manuscript in preparation.)

Please cite the method and this package if you use PALM-mbQTL in your work.
