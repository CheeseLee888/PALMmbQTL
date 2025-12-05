#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# options(stringsAsFactors=F, scipen = 999)
options(stringsAsFactors = F)
library(SAIGEQTL)


BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
print(sessionInfo())

option_list <- list(
  make_option("--inFile",
    type = "character", default = "",
    help = "Path to geno file. Support many formats including bed/bim/fam, bgen, gds and txt. If using bed/bim/fam, only go with prefix; Otherwise full file name."
  ),
  make_option("--GMMATmodelFile",
    type = "character", default = "",
    help = "Path to the input file containing the glmm model, which is output from previous step. Will be used by load()"
  ),
  make_option("--SAIGEOutputFile",
    type = "character", default = "",
    help = "Path to the output file containing assoc test results"
  )
)


parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


SPAGMMATtest(
  inFile = opt$inFile,
  GMMATmodelFile = opt$GMMATmodelFile,
  SAIGEOutputFile = opt$SAIGEOutputFile
)
