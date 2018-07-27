env_shiny <- new.env()

env_shiny$peak_merge2 <- function(peak, n = 250) {
  peak <- GenomicRanges::makeGRangesFromDataFrame(peak)
  GenomicRanges::end(peak) <- GenomicRanges::end(peak) + n
  peak <- GenomicRanges::reduce(peak)
  GenomicRanges::end(peak) <- GenomicRanges::end(peak) - n
  peak
}
