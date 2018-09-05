#' transform one hg coordinates to another
#'
#'
#' Convert coordinates between different genomes when necessary.
#'
#'
#' @param peak A \code{GRange} object. The genome, the coordinates of which need
#'   to be coverted.
#' @param hg_from The genome are coverting from. This parameter can be "hg18",
#'   "hg19" or "hg38", etc.
#' @param hg_to Which type the genome is converting to. This parameter can be "hg18",
#'   "hg19" or "hg38", etc. NOTICE \code{hg_from} and \code{hg_to} should be
#'   different from each other.
#'
#' @return A GRanges object.
#'
#' @examples
#'
#' \dontrun{
#' peak_df <- data.frame(chrom = c("chr7", "chr7", "chr7"),
#'                       chromStart = c(128043908, 128045075, 128046242),
#'                       chromEnd = c(128045074, 128046241, 128047408),
#'                       stringsAsFactors = FALSE)
#' peak <-  as(peak_df, "GRanges")
#'
#' trans_cor(peak, "hg19", "hg38")
#'
#' }
#'
#' @importFrom stringr str_to_title
#' @importFrom R.utils gunzip
#' @importFrom rtracklayer import.chain
#' @importFrom rtracklayer liftOver
#' @importFrom utils download.file
#'
#' @export

trans_cor <- function(peak, hg_from, hg_to) {
  file_name <- sprintf("%sTo%s.over.chain", hg_from, str_to_title(hg_to))
  path <- file.path(system.file(package="primirTSS", "extdata"), file_name)

  if (!file.exists(path)) {
    gz_name <- file.path(system.file(package="primirTSS", "extdata"),
                         paste0(file_name, ".gz"))
    url <- sprintf(
      "http://hgdownload.cse.ucsc.edu/goldenPath/%s/liftOver/%s.gz",
      hg_from, file_name)
    .trans_hg_download(url, gz_name)
    gunzip(gz_name)
  }
  path <- system.file(package="primirTSS", "extdata", file_name)
  ch = import.chain(path)
  suppressWarnings(unlist(liftOver(peak, ch)))
}


.trans_hg_download <- function(url, gz_name, download.file, N.TRIES=3L) {
  N.TRIES <- as.integer(N.TRIES)
  stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

  while (N.TRIES > 0L) {
    result <- tryCatch(utils::download.file(url, gz_name), error=identity)
    if (!inherits(result, "error"))
      break
    N.TRIES <- N.TRIES - 1L
  }

  if (N.TRIES == 0L) {
    stop("'trans_cor()' failed:",
         "\n  URL: ", url,
         "\n  error: ", conditionMessage(result))
  }

  result
}





