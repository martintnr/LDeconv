#' Compute deconvolved test statistics from LD inverse matrices and formatted summary statistics
#'
#' Applies the deconvolution step to a single chromosome-block unit by combining
#' the central block with its immediate neighbors, loading the corresponding LD
#' inverse matrix, and computing deconvolved test statistics.
#'
#' @param X Integer index of the chromosome-block unit to process.
#'
#' @param DATA A list of formatted data frames of summary statistics, as produced
#'   by `sumstats_deconvolution()`.
#'
#' @param ld_inverse_dir Path to the directory containing the LD inverse matrices.
#' These files can be downloaded from Zenodo at `zenodo/abcd.fr`.
#'
#' @returns A data frame for the processed block containing, for each variant,
#'   the variant identifier, the original test statistic, and the deconvolved
#'   test statistic. Returns `NULL` if the block contains no observed test
#'   statistics or if a matrix dimension mismatch is encountered.
#'
deconvolution_computation <- function(X, DATA, ld_inverse_dir){

  LDinv <- NULL
  t_py  <- NULL
  z_py  <- NULL

  on.exit({
    rm(LDinv, t_py, z_py)
    gc(full = TRUE)
    py_gc$collect()
    rm(list = ls())
    gc()
  }, add = TRUE)
 # Release the LD inverse matrix on from memory to limit RAM usage, on exit.



  Central <- DATA[[X]]
  # Extract the central chromosome-block unit to be processed.

  if(sum(is.na(Central$tstat))==nrow(Central)){return(NULL)}
  # Return `NULL` if the central block contains only missing test statistics.

  if(X>1){Before <- DATA[[X-1]]}else{Before <- NULL}
  if(X<length(DATA)){After <- DATA[[X+1]]}else{After <- NULL}

  # Retrieve the immediately preceding and following blocks, when available.

  if(!is.null(Before)){
    if(Before$chromosome[1] != Central$chromosome[1]){
      Before <- NULL}}
  if(!is.null(After)){
    if(After$chromosome[1] != Central$chromosome[1]){
      After <- NULL}}

  # Discard neighboring blocks if they belong to a different chromosome.

  Window <- rbind(Before, Central, After)
  Window$tstat[is.na(Window$tstat)] <- 0

  # Construct the full analysis window and replace missing test statistics with zero so that they do not contribute to the deconvolution.
  PATH <- paste0(ld_inverse_dir, "/chr", Central$chromosome[1],"_LD_inv_BD_", Central$block[1], ".npz")
  LDinv = scipy_sparse$load_npz(PATH)
  # Load the LD inverse matrix corresponding to the current central block.


  if(dim(LDinv)[1] != length(Window$tstat) | dim(LDinv)[2] != length(Window$tstat)){
    message(paste0("Different dimensions between inverse matrix and tstat vector for chromosome-block pair number ", X, " ; this should not happen"))
    return(NULL)
  }

  t_py <- r_to_py(as.numeric(Window$tstat))
  z_py <- LDinv$dot(t_py)
  Window$tstat_LDeconv <- as.numeric(py_to_r(z_py))
  # Compute the deconvolved test statistics by matrix-vector multiplication, through python


  Central <- merge(Central, Window[,c("variant", "tstat_LDeconv")], by = c("variant"), sort = F)
  # Retain only deconvolved test statistics corresponding to variants in the central block, which are the target output of the computation.


  message(paste0("Chromosome-block pair number ", X, " processed, out of ", length(DATA) , " total"))
  return(Central)

}
