env_shiny <- new.env()

env_shiny$peak_merge2 <- function(peak, n = 250) {
  peak <- GenomicRanges::makeGRangesFromDataFrame(peak)
  GenomicRanges::end(peak) <- GenomicRanges::end(peak) + n
  peak <- GenomicRanges::reduce(peak)
  GenomicRanges::end(peak) <- GenomicRanges::end(peak) - n
  peak
}

env_shiny$find_tss2 <- function(bed_merged, expressed_mir = "all",
                                flanking_num = 1000, threshold = 0.7,
                                ignore_DHS_check = TRUE,
                                DHS, allmirdhs_byforce = TRUE,
                                expressed_gene = "all", allmirgene_byforce = TRUE,
                                seek_tf = FALSE, tf_n = 1000, min.score = 0.8){
  DHS <- GenomicRanges::makeGRangesFromDataFrame(DHS)
  find_tss(bed_merged, expressed_mir, flanking_num,
           threshold, ignore_DHS_check, DHS = DHS, allmirdhs_byforce,
           expressed_gene, allmirgene_byforce, seek_tf,
           tf_n, min.score)
}

env_shiny$plot_primiRNA2 <- function(expressed_mir, bed_merged,
                                     flanking_num = 1000, threshold = 0.7,
                                     ignore_DHS_check = TRUE, DHS,
                                     allmirdhs_byforce = TRUE,
                                     expressed_gene = "all",
                                     allmirgene_byforce = TRUE){
  DHS <- GenomicRanges::makeGRangesFromDataFrame(DHS)
  plot_primiRNA(expressed_mir, bed_merged, flanking_num, threshold,
                ignore_DHS_check, DHS, allmirdhs_byforce,
                expressed_gene = "all", allmirgene_byforce = TRUE)
}

