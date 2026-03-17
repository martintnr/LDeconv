
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDeconv: obtaining LD-adjusted summary statistics through linkage disequilibrium deconvolution

<!-- badges: start -->

<!-- badges: end -->

## Principle of LDeconv

![LDeconv](Github_Fig.png) LDeconv takes as input GWAS summary
statistics of a trait and outputs deconvoluted LD-adjusted summary
statistics. Further information can be found in our
[pre-print](NotOnlineYet), coming soon.

## Installation

To install the current version of LDeconv, which should take a few
minutes:

``` r
if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

if(!"LDeconv" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/LDeconv") }
```

PRISM has been tested on Linux: Ubuntu 22.04.4 LTS and macOS: Tahoe 26.3

## Toy example

This example illustrates how to compute LD-adjusted summary statistics
for variants in the APOE region using UK Biobank round 2 triglyceride
summary statistics.

The script below creates a working directory and downloads from Zenodo
the LDeconv variant index together with the LD inverse matrices required
for the analysis. If computational resources permit, the `parallel` and
`nbcores` arguments can be used to enable parallel execution.

``` r
library(R.utils)
library(data.table)
library(reticulate)
#library(parallel)
py_require("scipy")
scipy_sparse <- import("scipy.sparse", convert = FALSE)
py_gc <- import("gc", convert = FALSE)


if(!grepl("LDeconv_ToyExample", getwd(), fixed = TRUE)){
  if(!file.exists("LDeconv_ToyExample/")){system("mkdir LDeconv_ToyExample")}
  setwd("LDeconv_ToyExample")}


ld_inverse_dir = "LDeconv_ToyExample"
index <- fread("LDeconv_ToyExample/Index_Variant_Block.csv")
ListofTraits <- unique(c(ParametersTable$X, ParametersTable$Y))

sumstats <- fread(LDeconv_ToyExample/Index_Variant_Block.csv)

Results <- LDeconv::sumstats_deconvolution(sumstats, ld_inverse_dir, index, parallel=F, nbcores = 1)
```

When run on a single CPU core, `sumstats_deconvolution()` is expected to
complete in a few minutes. The `Results` data frame contains the
original and LD-adjusted summary statistics for each processed variant.

To visualize a result manhattan plot[example
graph](https://github.com/martintnr/LDeconv/NotMadeYet):

``` r
library(ggplot2)

Graph <- ggplot()

print(Graph)
```

## Using PRISM

To use PRISM on real data, please use this
[vignette](https://github.com/martintnr/LDeconv/blob/main/vignettes/NotMadeYet).
