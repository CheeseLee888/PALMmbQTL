# PALM-mbQTL: Docker quick start

PALM-mbQTL implements a Poisson mixed-effects model for microbiome QTL
(mbQTL) analysis. For each microbial feature (e.g. genus-level
abundance), the pipeline

1. builds a genetic relationship matrix (GRM) from SNP genotypes;
2. merges microbiome abundance with sample-level covariates;
3. fits a null Poisson GLMM with random genetic effect; and
4. performs single-variant score tests across the genome.

This repository does **not** ship the source code of PALM-mbQTL.
Instead, it provides:

- a pre-built Docker image on Docker Hub, and
- a minimal example folder with toy data and shell wrappers,

so that users can quickly run the pipeline end-to-end without
installing R packages or compiling anything locally.

---

## Repository contents

- `quick_user_guide.pdf`  
  A short PDF describing the Docker-based installation and the example
  workflows.

- `example/`  
  Toy mbQTL dataset and shell scripts:
  - `input/` – small PLINK genotype, microbiome abundance, and
    covariate tables.
  - `output/` – will be created when the examples are run.
  - `quick.sh` – mbQTL pipeline (GRM + null GLMM + score test).
  - `quick_palm.sh` – original PALM framework (no genotype, only
    abundance vs covariates).

For the full statistical details of the PALM framework and all PALM
options, please refer to the separate PALM manual PDF
(e.g. `PALM_manual.pdf`) distributed with the software.

---

## Requirements

- A working Docker installation
  - Docker Desktop on macOS/Windows, or Docker Engine on Linux.
- Permission to run `docker` from the command line.

No local R installation is required; all dependencies are bundled
inside the Docker image.

---

## Pull the Docker image

Pull the pre-built PALM-mbQTL image from Docker Hub:

```bash
docker pull --platform=linux/amd64 cheeselee/palmmbqtl:latest
```

On a standard x86_64 Linux machine you may omit the `--platform`
option, but it is safe to keep it (and recommended for macOS with
Apple Silicon).

---

## Run the mbQTL example

From the `example/` directory:

```bash
cd /path/to/example

docker run --rm --platform=linux/amd64   -v "$PWD":/work   cheeselee/palmmbqtl:latest   bash -lc "cd /work && bash quick.sh"
```

This command:

- mounts the current directory into the container at `/work`,
- runs `quick.sh` inside the container, and
- writes all outputs to the local `output/` folder.

The example is intentionally tiny (a few samples and SNPs) so that it
finishes in a few seconds.

---

## Run the PALM example (no genotype)

Besides the Poisson mixed-effects mbQTL pipeline, the Docker image also
includes a built-in **PALM mode** for testing associations between
microbial features and sample-level covariates *without* using genotype
data.

From the same `example/` directory:

```bash
cd /path/to/example

docker run --rm --platform=linux/amd64   -v "$PWD":/work   cheeselee/palmmbqtl:latest   bash -lc "cd /work && bash quick_palm.sh"
```

Inputs are the abundance and covariate tables under `input/`, and
outputs are written to `output/` with the prefix specified in
`quick_palm.sh`.

The detailed meaning of PALM arguments and all advanced options
(depth/prevalence filtering, multiple-testing correction, meta-analysis,
compositional correction, etc.) are documented in the PALM manual PDF,
not repeated here.

---

## Using your own data

The shell wrappers are designed so that you only need to:

1. **Prepare input files**

   In `example/input/`, replace the toy data with your own:

   - `geno.bed`, `geno.bim`, `geno.fam` – PLINK genotype files
     (for mbQTL analysis).
   - `abd.tsv` – microbiome abundance table, with one row per sample
     and one column per microbial feature.
   - `cov.tsv` – covariate table with a shared sample ID column.

2. **Edit the parameter block**

   At the top of `quick.sh` you will find lines like:

   ```bash
   inputFolder=input
   outputFolder=output

   genoFile=geno
   abdFile=abd.tsv
   covFile=cov.tsv
   sampleIDColinabdFile=sample_id
   sampleIDColincovFile=sample_id

   phenoCol=g_Blautia
   offsetCol=SeqDepth
   covarColList=age,sex
   ```

   Adapt these to your filenames and column names as needed
   (e.g. change `phenoCol`, `covarColList`, etc.).

   Similarly, at the top of `quick_palm.sh` you can set:

   ```bash
   inputFolder=input
   outputFolder=output

   abdFile=abd.tsv
   covFile=cov.tsv
   sampleIDCol=sample_id

   covariateInterest=sex
   covariateAdjust="age,PC1,PC2,batch"

   outPrefix=PALM_results
   ```

3. **Re-run the Docker commands**

   Run the same `docker run ... quick.sh` or `quick_palm.sh` commands as
   above; the pipeline will now use your real data.

For more complex deployment (e.g. splitting the steps and submitting
them as separate jobs on an HPC scheduler, or integrating with existing
Docker/Apptainer/Singularity `.sif` workflows), please contact the
author.

---

## Citation and contact

If you use PALM-mbQTL or the PALM mode in your work, please cite the
corresponding method papers (mbQTL / PALM) once they are available.

For questions, bug reports, or help with custom deployments, please
contact:

- **Maintainer:** Peter Li
- **Email:** cheese.lee.888@gmail.com  
- **Docker Hub:** `cheeselee/palmmbqtl`
