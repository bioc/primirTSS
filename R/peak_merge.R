#' Merge adjacent peaks within H3K4me3 or Pol II data.
#'
#' Merge the adjacent segments provided as \code{GRange} object from original
#' data. This function will merge adjacent peaks the distance between which is
#' less than n base pairs apart and then return the merged segments.
#'
#' @param peak A \code{GRange} object. The peaks to be merged from one certain
#'   ChIP-seq data, such as H3K4me3 data or Pol II data.
#' @param n A number. \code{n} stipulates the distance(bp, base pair) between
#'   two seperate peaks within which they should be merged.
#'
#' @return A GRanges object. The merged peaks for the following analysis to
#'   search for TSSs.
#'
#'
#' @examples
#' peak_df <- data.frame(chrom = c("chr1", "chr2", "chr1"),
#'                       chromStart = c(450, 460, 680),
#'                       chromEnd = c(470, 480, 710),
#'                       stringsAsFactors = FALSE)
#' peak <-  as(peak_df, "GRanges")
#'
#' peak_merge(peak, n =250)
#'
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer chrom
#' @importFrom BiocGenerics strand
#' @importFrom Gviz seqnames
#' @importFrom Gviz symbol
#' @export

peak_merge <- function(peak, n = 250) {
  end(peak) <- end(peak) + n
  peak <- reduce(peak)
  end(peak) <- end(peak) - n
  peak
}
