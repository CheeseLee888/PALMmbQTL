---
title: Installation
nav_order: 3
---

<!--
This page explains how to install PALM-mbQTL (SAIGEQTL) and set up
its software environment. Adjust URLs and commands to match your repo.
-->

# Installation

PALM-mbQTL is implemented as an R package (internal name: `SAIGEQTL`)
with additional helper scripts in the `extdata/` directory.

You can run PALM-mbQTL using:

- a local pixi environment (recommended),
- a Docker image, or
- an Apptainer/Singularity image built from the Docker image.

## Requirements

- Linux or macOS (for full support)
- R version 3.5.0 or higher
- C++14 compiler with OpenMP support
- GNU make

The Docker/Apptainer images bundle all required R dependencies.

## Option 1: pixi-based installation (local development)

```bash
# 1. Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"

# 2. Clone the PALMmbQTL repository
git clone https://github.com/USERNAME/PALMmbQTL.git
cd PALMmbQTL

# 3. Create the environment and install R dependencies
pixi install

# 4. Install the SAIGEQTL package from the current directory
pixi run R CMD INSTALL .
```

After this, you can run R scripts inside the pixi environment using:

```bash
pixi run --manifest-path=./pixi.toml Rscript script.R
```

## Option 2: Docker image

You can build a local Docker image from the repository:

```bash
cd /path/to/PALMmbQTL
docker build -t palmmbqtl:latest .
```

Test that the wrapper scripts and package are available:

```bash
docker run --rm palmmbqtl:latest step1_fitNULL.R --help
```

Later, you can mount your own data into the container, for example:

```bash
docker run --rm -v /path/to/mydata:/data palmmbqtl:latest   step1_fitNULL.R --phenoFile=/data/pheno.txt ...
```

## Option 3: Apptainer / Singularity (HPC clusters)

Many HPC clusters do not allow Docker but support Apptainer/Singularity.
You can build a `.sif` image from the Docker image:

```bash
# On a machine where Docker is available and the palmmbqtl:latest
# image has been built or pulled:
apptainer build palmmbqtl_latest.sif docker-daemon://palmmbqtl:latest
```

Then copy `palmmbqtl_latest.sif` to your cluster and run:

```bash
apptainer exec palmmbqtl_latest.sif step1_fitNULL.R --help
```

To access external data, bind a directory:

```bash
apptainer exec --bind /project/mydata:/data palmmbqtl_latest.sif   step1_fitNULL.R --phenoFile=/data/pheno.txt ...
```

Once the environment is installed, proceed to the
[Quick Start guide](quickstart.md).
