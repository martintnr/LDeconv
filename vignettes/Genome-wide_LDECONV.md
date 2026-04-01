If you are reading this vignette, you likely want to apply LDeconv to
your own summary statistics. This tutorial explains how to do that step
by step.

Conceptually, the workflow is fairly simple. With the functions
provided, the core operation is essentially matrix multiplication. In
practice, the main challenge is downloading, organizing, and using the
genome-wide inverse LD matrices required for the analysis.

In this vignette, we illustrate the procedure by deconvoluting summary
statistics for chromosome 19 of triglycerides from UK Biobank GWAS round
2. The goal is to provide a reproducible example that you can easily
adapt to your own data.

1.  Standard setup: load the packages required to run LDeconv and to
    read the data.

``` r
library(LDeconv)
library(R.utils)
library(data.table)
library(reticulate)
#library(parallel)
library(LDeconv)

py_require("scipy")
scipy_sparse <- import("scipy.sparse", convert = FALSE)
py_gc <- import("gc", convert = FALSE)
# Installing python dependencies
```

1.  File download

Here, I create a working directory and the required subfolders, then
download the necessary files:

- the variant-LD index,
- the summary statistics for triglycerides from UK Biobank GWAS round 2
  (our example),
- the inverse LD matrices.

In this example, the inverse LD matrices are downloaded only for
chromosome 19, but the code can easily be adapted to download all
chromosomes.

``` r


if(!grepl("LDeconv_GenomeWide", getwd(), fixed = TRUE)){
  if(!file.exists("LDeconv_GenomeWide/")){system("mkdir LDeconv_GenomeWide")}
  setwd("LDeconv_GenomeWide")}

system("mkdir Sumstats")
system("wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz -O Sumstats/30870_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")

system("mkdir Index")
system("wget https://zenodo.org/records/19183010/files/Index_LDeconv.csv -O Index/Index_LDeconv.csv")


system("mkdir Inverse_Matrices")
#for(chromosome in c(1:22)){ 
for(chromosome in c(19)){
system(paste0("wget https://zenodo.org/records/19183010/files/chr",chromosome,"_npz.zip -O Inverse_Matrices/chr",chromosome,"_npz.zip"))
system(paste0("unzip Inverse_Matrices/chr",chromosome,"_npz.zip -d Inverse_Matrices"))
system(paste0("rm Inverse_Matrices/chr",chromosome,"_npz.zip"))
}

system("mkdir Results")
```

1.  Code Execution

We can now load the data, run the deconvolution, and save the results.

For this example, we first restrict the summary statistics to variants
on chromosome 19 and retain only the two columns required for
deconvolution.

If one or more required inverse LD matrices are missing, the analysis
will not complete. In contrast, missing variants in your summary
statistics will not prevent the function from running. Variants that are
absent from the index will simply be filtered out internally.

Keep in mind that memory usage can be substantial, especially when using
parallelization.

``` r


ld_inverse_dir = "Inverse_Matrices"
index <- fread("Index/Index_LDeconv.csv")
sumstats <- fread("Sumstats/30870_irnt.gwas.imputed_v3.both_sexes.tsv.bgz")

sumstats <- sumstats[grepl("^19:", variant), c("variant", "tstat")]

Results <- LDeconv::sumstats_deconvolution(sumstats, ld_inverse_dir, index, parallel=F, nbcores = 1)

write.csv(Results, "Results/VignetteExample.csv", row.names = F)
```
