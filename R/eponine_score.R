#' @importFrom dplyr mutate
#' @importFrom dplyr data_frame
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' @importFrom tidyr unite
#' @importFrom tidyr separate
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom Biostrings reverseComplement
#' @importFrom purrr pmap
#' @importFrom stats runif

cor_tdga <- function(chrom, start, end, strand) {
  a <- Hsapiens[[chrom]]
  if (strand == "+") {
    a[start:end] %>%
      as.character()
  } else {
    a[start:end] %>%
      reverseComplement() %>%
      as.character()
  }
}


require_fa <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                       strand, peak_p1, peak_p2, histone_p1_flank, histone_p2_flank) {

  mir_flank <- data_frame(mir_name = mir_name, chrom = chrom, stem_loop_p1 = stem_loop_p1,
                          stem_loop_p2 = stem_loop_p2, strand = strand,
                          peak_p1 = peak_p1, peak_p2 = peak_p2,
                          histone_p1_flank = histone_p1_flank,
                          histone_p2_flank = histone_p2_flank)
  tcga <- mir_flank %>%
    select(chrom, start = histone_p1_flank, end = histone_p2_flank, strand) %>%
    pmap(cor_tdga) %>%
    unlist()

  line1 <- mir_flank %>%
    mutate(arrow = ">") %>%
    unite(line1,
          mir_name, chrom, stem_loop_p1, stem_loop_p2, strand,
          peak_p1, peak_p2,
          histone_p1_flank, histone_p2_flank, sep = "_") %>%
    unite(line1, arrow, line1, sep = "")

  paste(line1$line1, tcga, sep = "\n")
}


eponine_score <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                          strand, peak_p1, peak_p2, flanking_num = 1000, threshold = 0.7) {

  mir_peaks <- data_frame(mir_name = mir_name, chrom = chrom, stem_loop_p1 = stem_loop_p1,
                          stem_loop_p2 = stem_loop_p2, strand = strand,
                          peak_p1 = peak_p1, peak_p2 = peak_p2)

  mir_flank <- mir_peaks %>%
    mutate(histone_p1_flank = peak_p1 - flanking_num,
           histone_p2_flank = peak_p2 + flanking_num)

  a <- require_fa(mir_flank$mir_name, mir_flank$chrom,
                  mir_flank$stem_loop_p1, mir_flank$stem_loop_p2,
                  mir_flank$strand, mir_flank$peak_p1, mir_flank$peak_p2,
                  mir_flank$histone_p1_flank, mir_flank$histone_p2_flank)
  tmp_path <- file.path(tempdir(), paste0(runif(1, 0, 100), ".fa"))
  writeLines(a, tmp_path)

  java_path <- system.file("extdata", "eponine-scan.jar", package = "primirTSS")
  cmd <- sprintf("java -jar %s -seq %s -threshold %s", java_path, tmp_path, threshold)
  aa <- system(cmd, intern = TRUE)

  file.remove(tmp_path)

  a_tmp <- aa[-(1:3)] %>%
    str_split("\t", simplify = TRUE)

  result <- data_frame(previous = a_tmp[, 1],
                       tss_p1 = a_tmp[, 4],
                       tss_p2 = a_tmp[, 5],
                       eponine_score = a_tmp[, 6]) %>%
    separate(previous,
             into = c("mir_name", "chrom", "stem_loop_p1", "stem_loop_p2", "strand",
                      "peak_p1", "peak_p2",
                      "histone_p1_flank", "histone_p2_flank"),
             sep = "_") %>%
    mutate(stem_loop_p1 = as.double(stem_loop_p1),
           stem_loop_p2 = as.double(stem_loop_p2),
           peak_p1 = as.double(peak_p1),
           peak_p2 = as.double(peak_p2),
           histone_p1_flank = as.double(histone_p1_flank),
           histone_p2_flank = as.double(histone_p2_flank),
           tss_p1 = as.double(tss_p1),
           tss_p2 = as.double(tss_p2),
           eponine_score = as.double(eponine_score)) %>%
    mutate(tss_p1 = histone_p1_flank + tss_p1 - 1,
           tss_p2 = histone_p1_flank + tss_p2 - 1) %>%
    filter(tss_p1 >= peak_p1 & tss_p2 <= peak_p2) %>%
    select(-(peak_p1:histone_p2_flank))

  fail_mir <- setdiff(unique(mir_peaks$mir_name), unique(result$mir_name))

  if (length(fail_mir) == 0) {
    fail <- "All miRNAs have eponine socres"
  } else {
    fail <- fail_mir
  }

  list(success = result, fail_eponine = fail)
}
