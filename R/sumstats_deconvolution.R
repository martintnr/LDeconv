#' Format summary statistics and run the deconvolution procedure
#'
#' Prepares input summary statistics for the deconvolution algorithm
#' and maps variants to the corresponding entries in the LD inverse matrices.
#'
#' @param sumstats A data frame containing exactly two columns:
#'   \describe{
#'     \item{variant}{Character vector of variant identifiers in the format
#'     `CHR:POS:A1:A2`, based on hg19 genomic coordinates.}
#'     \item{tstat}{Numeric vector of test statistics to be deconvolved.}
#'   }
#'
#' @param ld_inverse_dir Path to the directory containing the LD inverse matrices.
#' These files can be downloaded from Zenodo at `zenodo/abcd.fr`.
#'
#' @param index A variant index based on hg19 coordinates, used to match the
#'   variants in `sumstats` to their corresponding positions in the LD inverse
#'   matrices. This file can be downloaded from Zenodo at `zenodo/abcd.fr`.
#'
#' @param parallel Indicating whether block-wise computation
#'   should be performed in parallel using `parallel::mclapply()`.
#'   Parallel computation might be RAM expensive.
#'
#' @param nbcores Number of CPU cores to use when `parallel = TRUE`.
#'
#' @returns A data frame containing, for each processed variant, the index,
#'   the original test statistic, and the deconvolved test statistic.
sumstats_deconvolution <- function(sumstats, ld_inverse_dir, index, parallel=F, nbcores = 1){

  DATA <- merge(index, sumstats, by = "variant",
                 all.x = TRUE, all.y = F, sort = F)
  DATA <- as.data.frame(DATA)
  if(sum(DATA$variant == index$variant)!=nrow(DATA)){message("PROBLEM: index and data mismatch")}

  # Construct the summary statistics matrix for deconvolution after matching input variants to the variant index.
  # Variants present in `sumstats` but absent from `index` are excluded.

  key <- interaction(DATA$chromosome, DATA$block, drop = TRUE)
  pairs <- unique(DATA[c("chromosome", "block")])
  pairs <- pairs[order(pairs$chromosome, pairs$block), ]
  pairs$hit <- tapply(!is.na(DATA$tstat), key, any)[
    interaction(pairs$chromosome, pairs$block, drop = TRUE)]
  pairs$hit[is.na(pairs$hit)] <- FALSE
  pairs$keep <- with(pairs, ave(hit, chromosome, FUN = function(z) {
    z | c(FALSE, z[-length(z)]) | c(z[-1], FALSE)
  }))
  DATA <- DATA[
    key %in% interaction(pairs$chromosome[pairs$keep],
                         pairs$block[pairs$keep], drop = TRUE),]


  # Retain only chromosome-block units containing at least one observed test statistic, together with their immediately adjacent blocks.

  DATA <- split(DATA, interaction(DATA$block, DATA$chromosome))
  DATA <- DATA[sapply(DATA, nrow) > 0]
  # Split the data by chromosome-block unit for block-wise computation.


  if(parallel == TRUE){
  DATA <- mclapply(X = c(1:length(DATA)), FUN = deconvolution_computation, DATA=DATA, ld_inverse_dir, mc.cores=nbcores)
  }else{ DATA <- lapply(X = c(1:length(DATA)),FUN = deconvolution_computation, DATA=DATA, ld_inverse_dir)
  }
  # Perform block-wise deconvolution, optionally in parallel across `nbcores` using `parallel::mclapply()`.

  DATA <- do.call("rbind", DATA)
  return(DATA)

}












