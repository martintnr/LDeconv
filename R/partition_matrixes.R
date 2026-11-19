#' Identify and refine linkage disequilibrium blocks for each chromosome
#'
#' Reads a sparse LD correlation matrix for each
#' chromosome. Variants connected by correlations whose absolute values exceed
#' the threshold form initial blocks. The function then
#' splits oversized blocks using adjacency-constrained
#' hierarchical clustering.
#'
#' @param input_dir Path to the directory containing the upper-triangular LD
#'   matrices in SciPy sparse `.npz` format.

#'
#' @param output_dir Path to directory in which to save
#'   the block labels.
#'
#' @param threshold Numeric cutoff for the absolute LD correlation between two
#'   variants.
#'
#' @param max_block_size Positive integer giving the target maximum number of
#'   variants per block and the bandwidth used for clustering oversized
#'   blocks.
#'
#' @param chromosomes Vector of chromosome numbers to be processed.
#'
#' @returns On successful completion, returns `NULL`. For each
#'   chromosome, writes `chr<N>_refined_block_labels.npy` to `output_dir`.
#' @export
sumstats_deconvolution <- function( input_dir = "/Folder/With/Sparse/LDMatrix",
                                    output_dir = "Output/Directory",
                                    threshold = 1e-2,
                                    max_block_size = 10000,
                                    chromosomes = c(22:1)){

  library(reticulate)
  library(adjclust)
  library(Matrix)
  library(igraph)

  py_require("scipy")
  scipy_sparse <- import("scipy.sparse")


  # Identify LD blocks from thresholded upper-triangular correlations
  identify_blocks <- function(LD, threshold = 1e-2) {
    entries <- Matrix::summary(LD)
    keep <- (entries$i < entries$j) & (entries$x > threshold)

    graph <- igraph::make_empty_graph(n = nrow(LD), directed = FALSE)
    if (any(keep)) {
      edges <- as.vector(rbind(entries$i[keep], entries$j[keep]))
      graph <- igraph::add_edges(graph, edges)
    }

    labels <- igraph::components(graph)$membership
    split(seq_len(nrow(LD)), labels)
  }


  # Split initial components into independent LD blocks of at most 1,000,000 SNPs
  refine_blocks <- function(LD, blocks, max_block_size = 1000000) {
    refined_labels <- rep(NA_integer_, nrow(LD))
    new_block_id <- 1L

    for (block_id in sort(as.integer(names(blocks)))) {
      indices <- as.integer(blocks[[as.character(block_id)]])

      for (start in seq.int(1L, length(indices), by = max_block_size)) {
        stop <- min(start + max_block_size - 1L, length(indices))
        refined_labels[indices[start:stop]] <- new_block_id
        new_block_id <- new_block_id + 1L
      }
    }

    refined_labels
  }


  # Merge interrupted runs and attach singletons to longer neighboring runs
  smooth_interruptions <- function(labels) {
    runs <- rle(labels)
    n <- length(runs$lengths)

    for (i in seq_len(n)) {
      if (i > 1 && i < n &&
          runs$values[i - 1] == runs$values[i + 1] &&
          runs$values[i] != runs$values[i - 1]) {
        runs$values[i] <- runs$values[i - 1]
      } else if (runs$lengths[i] == 1) {
        left_len <- if (i > 1) runs$lengths[i - 1] else -Inf
        right_len <- if (i < n) runs$lengths[i + 1] else -Inf

        if (left_len >= right_len && left_len > 1) {
          runs$values[i] <- runs$values[i - 1]
        } else if (right_len > left_len && right_len > 1) {
          runs$values[i] <- runs$values[i + 1]
        }
      }
    }

    smoothed <- inverse.rle(runs)
    match(smoothed, unique(smoothed))
  }


  # Represent an upper-triangular block as a symmetric sparse matrix
  upper_to_symmetric <- function(LD) {
    entries <- Matrix::summary(LD)
    keep <- entries$i <= entries$j

    Matrix::sparseMatrix(
      i = entries$i[keep],
      j = entries$j[keep],
      x = entries$x[keep],
      dims = dim(LD),
      symmetric = TRUE
    )
  }


  # Split large blocks with adjacency-constrained clustering
  cluster_large_blocks <- function(block_id, LD, labels, max_block_size = 10000) {
    indices <- which(labels == block_id)
    block <- LD[indices, indices, drop = FALSE]
    block <- upper_to_symmetric(block)

    tree <- adjClust(
      block,
      type = "similarity",
      h = max_block_size,
      strictCheck = FALSE,
      nthreads = 5L
    )

    n_clusters <- 2L
    repeat {
      cluster_labels <- cutree_chac(tree, k = n_clusters)
      if (max(table(cluster_labels)) <= max_block_size) {
        break
      }
      n_clusters <- n_clusters + 1L
    }

    cluster_labels + (block_id * 1000)
  }


  # Main execution starts here


  process_chromosomes <- function(
      input_dir,
      output_dir,
      threshold = threshold,
      max_block_size = max_block_size,
      chromosomes = chromosomes
  ) {
    for (chromosome in chromosomes) {
      cat(sprintf("\nProcessing chromosome %d...\n", chromosome))

      ld_path <- file.path(
        input_dir,
        sprintf("chr%d_sparse_matrix_upper.npz", chromosome)
      )
      LD <- abs(scipy_sparse$load_npz(ld_path))

      blocks <- identify_blocks(LD, threshold = threshold)
      refined_labels <- refine_blocks(LD, blocks, max_block_size = 1000000)
      refined_labels <- smooth_interruptions(refined_labels)

      block_sizes <- table(refined_labels)
      large_blocks <- as.integer(names(block_sizes)[block_sizes > max_block_size])

      if (length(large_blocks) > 0) {
        split_labels <- lapply(
          large_blocks,
          cluster_large_blocks,
          LD = LD,
          labels = refined_labels,
          max_block_size = max_block_size
        )

        for (i in seq_along(large_blocks)) {
          refined_labels[refined_labels == large_blocks[i]] <- split_labels[[i]]
        }

        refined_labels <- match(refined_labels, unique(refined_labels))
      }

      RcppCNPy::npySave(
        file.path(output_dir, sprintf("chr%d_refined_block_labels.npy", chromosome)),
        refined_labels - 1L,
        mode = "integer"
      )

      gc()
    }
  }


  process_chromosomes(
    input_dir = input_dir,
    output_dir = output_dir,
    threshold = threshold ,
    max_block_size = max_block_size,
    chromosomes = chromosomes)
}
