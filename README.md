
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LDeconv: obtaining LD-adjusted summary statistics through linkage disequilibrium deconvolution

<!-- badges: start -->

<!-- badges: end -->

## Principle of LDeconv

<!--   ![LDeconv](Github_Fig.png) -->

LDeconv takes as input GWAS summary statistics of a trait and outputs
deconvoluted LD-adjusted summary statistics. Further information can be
found in our [pre-print](NotOnlineYet), coming soon.

## Installation

To install the current version of LDeconv, which should take a few
minutes:

``` r
if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

if(!"LDeconv" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/LDeconv") }
```

LDeconv has been tested on Linux: Ubuntu 25.10 and macOS: Tahoe 26.3

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
library(LDeconv)

py_require("scipy")
scipy_sparse <- import("scipy.sparse", convert = FALSE)
py_gc <- import("gc", convert = FALSE)
# Installing python dependencies


if(!grepl("LDeconv_ToyExample", getwd(), fixed = TRUE)){
  if(!file.exists("LDeconv_ToyExample/")){system("mkdir LDeconv_ToyExample")}
  setwd("LDeconv_ToyExample")}
system("wget https://zenodo.org/records/19071985/files/chr19_LD_inv_BD_86.npz")
system("wget https://zenodo.org/records/19071985/files/Index_LDeconv.csv")
system("wget https://zenodo.org/records/19071985/files/ToyExample_Triglycerides.csv")
# Creates a new directory and download the necessary files from Zenodo, <200 MB


ld_inverse_dir = getwd()
index <- fread("Index_LDeconv.csv")
sumstats <- fread("ToyExample_Triglycerides.csv")

Results <- LDeconv::sumstats_deconvolution(sumstats, ld_inverse_dir, index, parallel=F, nbcores = 1)
```

When run on a single CPU core, `sumstats_deconvolution()` is expected to
complete in less than a minute. The `Results` data frame contains the
original and LD-adjusted summary statistics for each processed variant.

To visualize a result [manhattan
plot](https://github.com/martintnr/LDeconv/blob/main/ManhattanPlot_ToyExample.png):

``` r
#install.packages("ggplot2")
library(ggplot2)

# Compute two-sided p-values from t statistics
Results$p_gwas <- 2 * pnorm(-abs(Results$tstat))
Results$p_ld   <- 2 * pnorm(-abs(Results$tstat_LDeconv))

# Avoid small p-value for the graph
Results$p_gwas <- pmax(Results$p_gwas, 1e-20)
Results$p_ld   <- pmax(Results$p_ld, 1e-20)

# Reshape to long format for ggplot
plot_Results <- rbind(
  data.frame(position = Results$position,
    chromosome = Results$chromosome,
    variant = Results$variant,
    method = "Original GWAS",
    p = Results$p_gwas),
  data.frame(position = Results$position,
    chromosome = Results$chromosome,
    variant = Results$variant,
    method = "LD-adjusted",
    p = Results$p_ld)
)

# Order legend
plot_Results$method <- factor(
  plot_Results$method,
  levels = c("Original GWAS", "LD-adjusted"))

# Manhattan y-axis in -log10
plot_Results$logp <- -log10(plot_Results$p)

# Genome-wide significance threshold
sig_line <- -log10(5e-8)

Graph <- ggplot(plot_Results, aes(x = position, y = logp, color = method, shape = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = sig_line, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Original GWAS" = "forestgreen",
    "LD-adjusted"   = "dodgerblue3")) +
  scale_shape_manual(values = c(
    "Original GWAS" = 16,   
    "LD-adjusted"   = 17)) +
  labs(x = paste0("Genomic position on chromosome ", unique(Results$chromosome)),
    y = expression(-log[10](p)), color = NULL, shape = NULL) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom",panel.grid.minor = element_blank())

print(Graph)
```

## Using LDeconv

To use LDeconv on real data, please use this
[vignette](https://github.com/martintnr/LDeconv/blob/main/vignettes/NotMadeYet)
(coming soon !)
