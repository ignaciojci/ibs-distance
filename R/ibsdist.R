#' Compute IBS distance matrix using bit-parallel C++ implementation
#'
#' @param M Integer matrix with samples in rows and markers in columns
#' @param block_size Block size for cache optimization
#' @param diagonal_trueIBS Logical, if TRUE uses diagonal = 1
#' @param min_sites Minimum number of sites to compute distance
#' @param n_threads Number of threads (0 = auto)
#'
#' @return Numeric matrix of pairwise IBS distances
#' @export
ibs.dist <- function(M, block_size = 256, diagonal_trueIBS = FALSE,
                       min_sites = 0, n_threads = 0) {
    fast_ibs(M, block_size, diagonal_trueIBS, min_sites, n_threads)
}
