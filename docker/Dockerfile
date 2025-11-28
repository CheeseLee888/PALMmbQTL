FROM ubuntu:20.04

LABEL maintainer="cheese.lee.888@gmail.com"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    apt-get install -y pigz && \
    apt-get -y install curl time && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")' && \
    pixi run Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("barkasn/fastSave")' && \
    pixi run Rscript -e 'install.packages(c("GMMAT"), repos="https://cloud.r-project.org")'

RUN pixi run Rscript - << 'EOF'
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(
  c("gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData"),
  ask        = FALSE,
  update     = FALSE,
  dependencies = TRUE
)

ip <- rownames(installed.packages())
needed <- c("gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData")
missing <- setdiff(needed, ip)

if (length(missing)) {
  stop("Missing Bioconductor packages after install: ",
       paste(missing, collapse = ", "))
}
EOF

# delete any precompiled files to ensure clean build
RUN rm -f src/*.so src/*.o || true

## 2. Install SAIGEQTL package (key)
RUN pixi run R CMD INSTALL .

# Check all key packages at once
RUN pixi run Rscript - << 'EOF'
ip <- rownames(installed.packages())

needed <- c(
  "lintools", "remotes", "fastSave", "GMMAT",
  "gdsfmt", "SNPRelate", "SeqArray", "GenomeInfoDbData",
  "SAIGEQTL"
)
missing <- setdiff(needed, ip)
if (length(missing)) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "))
}

for (pkg in c("SAIGEQTL", "GMMAT", "SNPRelate", "SeqArray")) {
  message("Testing library(", pkg, ") ...")
  library(pkg, character.only = TRUE)
}
message("All key packages installed and loadable.")
EOF

# Create Singularity-compatible wrapper scripts
RUN for script in step0_generateGRM.R step0_mergePheno.R step1_fitNULL.R step2_scoreTest.R; do \
      echo '#!/bin/bash' > /usr/local/bin/$script && \
      echo 'cd /app' >> /usr/local/bin/$script && \
      echo 'exec pixi run Rscript /app/extdata/'$script' "$@"' >> /usr/local/bin/$script && \
      chmod +x /usr/local/bin/$script; \
    done


RUN mv extdata/input/ /usr/local/bin/

RUN rm -r src/ || true
RUN rm -r R/   || true

