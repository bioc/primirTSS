#' Integrate H3K4me3 data and Pol II data.
#'
#' Integrate peaks from H3K4me3 and Pol II data. To conduct the overlappd ranges
#' for the further analysis by imposing H3K4me3 peaks on Pol II peaks, if both
#' of these two different kinds of ChIP-seq data are available.
#'
#' @param peak1 H3K4me3 peaks. Merged peak data as \code{GRange} object by
#'   fucntion {peak_merge}
#' @param peak2 Pol II peaks. Merged peak data as \code{GRange} object by
#'   fucntion {peak_merge}
#'
#' @section Detail: Peak1 and peak2 are signals seperately from the ChIP-seq
#'   data of H3K4me3 and Pol II data that to be integrated. The data is
#'   \code{GRange} object containing three columns \code{Chrom}, \code{Ranges},
#'   \code{Strand}. And the order of these two kinds of data when input as
#'   \code{peak1} and \code{peak2} can be swapped.
#'
#' @seealso {peak_merge}
#'
#' @return A GRanges object. The joined peaks for the following analysis to
#'   search for TSSs.
#'
#' @examples
#' peak_df1 <- data.frame(chrom = c("chr1", "chr1", "chr1", "chr2"),
#'                        start = c(100, 460, 600, 70),
#'                        end = c(200, 500, 630, 100),
#'                        stringsAsFactors = FALSE)
#' peak1 <-  as(peak_df1, "GRanges")
#'
#' peak_df2 <- data.frame(chrom = c("chr1", "chr1", "chr1", "chr2"),
#'                        start = c(160, 470, 640, 71),
#'                        end = c(210, 480, 700, 90),
#'                        stringsAsFactors = FALSE)
#' peak2 <-  as(peak_df2, "GRanges")
#'
#' peak_join(peak1, peak2)
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges intersect
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export

peak_join <- function(peak1, peak2) {
  ov <- findOverlaps(peak1, peak2)
  intersect(peak1[queryHits(ov)], peak2[subjectHits(ov)])
}
