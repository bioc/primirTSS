#' @importFrom phastCons100way.UCSC.hg38 phastCons100way.UCSC.hg38
#' @importFrom GenomicScores gscores
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup


phast_score <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                        strand, tss_p1, tss_p2, eponine_score) {

  mir_eponine_score <- data_frame(mir_name = mir_name, chrom = chrom,
                                  stem_loop_p1 = stem_loop_p1,
                                  stem_loop_p2 = stem_loop_p2,
                                  strand = strand, tss_p1 = tss_p1,
                                  tss_p2 = tss_p2,
                                  eponine_score = eponine_score)

  phast <- phastCons100way.UCSC.hg38

  score_s <- function(loci) {
    gscores(phast, GRanges(loci))$default
  }

  mir_eponine_score %>%
    mutate(loci = paste(chrom, paste(tss_p1, tss_p2, sep = "-"), sep = ":"),
           phast_score = score_s(loci)) %>%
    group_by(mir_name) %>%
    mutate(eponine_rank = rank(eponine_score),
           phast_rank = rank(phast_score),
           e_p_rank = eponine_rank + phast_rank) %>%
    filter(e_p_rank == max(e_p_rank)) %>%
    ungroup() %>%
    select(-(loci:e_p_rank), -eponine_score)
}

phast_score_plot <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                             strand, tss_p1, tss_p2, eponine_score) {

  mir_eponine_score <- data_frame(mir_name = mir_name, chrom = chrom,
                                  stem_loop_p1 = stem_loop_p1,
                                  stem_loop_p2 = stem_loop_p2,
                                  strand = strand, tss_p1 = tss_p1,
                                  tss_p2 = tss_p2,
                                  eponine_score = eponine_score)

  phast <- phastCons100way.UCSC.hg38

  score_s <- function(loci) {
    gscores(phast, GRanges(loci))$default
  }

  mir_eponine_score %>%
    mutate(loci = paste(chrom, paste(tss_p1, tss_p2, sep = "-"), sep = ":"),
           phast_score = score_s(loci))
}



#' @importFrom IRanges IRanges
#' @importFrom IRanges countOverlaps
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#'
check_DHS_s <- function(chrom, strand, tss_p1, tss_p2, DHS) {

  DHS <- as.data.frame(DHS)
  DHS <- DHS %>%
    mutate(chrom = as.character(seqnames),
           dhs_p1 = as.double(start),
           dhs_p2 = as.double(end)) %>%
    select(chrom, dhs_p1, dhs_p2)

  if (strand == "+") {
    DHS_d <- DHS %>%
      mutate(dhs_p1_down = dhs_p1,
             dhs_p2_down = dhs_p2 + 1000)
  } else {
    DHS_d <- DHS %>%
      mutate(dhs_p1_down = dhs_p1 - 1000,
             dhs_p2_down = dhs_p2)
  }


  DHS_dc <- DHS_d[DHS_d$chrom == chrom, ]

  ir1 <- IRanges(start = tss_p1, end = tss_p2)
  ir2 <- IRanges(start = DHS_dc$dhs_p1_down, end = DHS_dc$dhs_p2_down)

  if (countOverlaps(ir1, ir2) != 0) {
    ov <- findOverlaps(ir1, ir2)
    tmp <- IRanges::intersect(ir1[queryHits(ov)], ir2[subjectHits(ov)])
    result <- paste(paste(start(tmp), end(tmp), sep = "_"), collapse = ":")
  } else {
    result <- "*"
  }

  result
}

check_DHS_df <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                         strand, tss_p1, tss_p2, eponine_score, DHS) {

  x <- data_frame(mir_name = mir_name, chrom = chrom,
                  stem_loop_p1 = stem_loop_p1,
                  stem_loop_p2 = stem_loop_p2,
                  strand = strand, tss_p1 = tss_p1,
                  tss_p2 = tss_p2,
                  eponine_score = eponine_score)
  candidate_tss <- x %>%
    select(chrom, strand, tss_p1, tss_p2) %>%
    pmap(check_DHS_s, DHS) %>%
    unlist()

  new_mir_tss_df1 <- x %>%
    bind_cols(can_tss = candidate_tss) %>%
    filter(can_tss != "*") %>%
    extand_df("can_tss", ":") %>%
    select(-tss_p1, -tss_p2) %>%
    separate(new_info, into = c("tss_p1", "tss_p2"), sep = "_") %>%
    mutate(tss_p1 = as.double(tss_p1),
           tss_p2 = as.double(tss_p2))

  new_mir_tss_df2 <- x %>%
    bind_cols(can_tss = candidate_tss)

  f2 <- function(x) {
    if (all(unique(x$can_tss) == "*")) {
      unique(x$mir_name)
    } else {
      NULL
    }
  }

  DHS_fail_mir <- lapply(split(new_mir_tss_df2, new_mir_tss_df2$mir_name),
                         f2) %>% unlist()

  list(success = new_mir_tss_df1, fail_DHS = DHS_fail_mir)
}


find_candidate_tss <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                               strand, tss_p1, tss_p2, eponine_score,
                               ignore_DHS_check = TRUE,
                               DHS, allmirdhs_byforce = TRUE) {

  mir_eponine_score <- data_frame(mir_name = mir_name, chrom = chrom,
                                  stem_loop_p1 = stem_loop_p1,
                                  stem_loop_p2 = stem_loop_p2,
                                  strand = strand, tss_p1 = tss_p1,
                                  tss_p2 = tss_p2,
                                  eponine_score = eponine_score)

  if (ignore_DHS_check == TRUE) {

    a <- phast_score(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                     strand, tss_p1, tss_p2, eponine_score)
    fail <- "Do not have a DHS check"
  } else {

    dhs_check <- check_DHS_df(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                              strand, tss_p1, tss_p2, eponine_score, DHS)

    mir_eponine_score_s <- dhs_check$success

    a <- phast_score(mir_eponine_score_s$mir_name, mir_eponine_score_s$chrom,
                     mir_eponine_score_s$stem_loop_p1,
                     mir_eponine_score_s$stem_loop_p2,
                     mir_eponine_score_s$strand,
                     mir_eponine_score_s$tss_p1, mir_eponine_score_s$tss_p2,
                     mir_eponine_score_s$eponine_score)

    if (is.null(dhs_check$fail_DHS)) {
      fail <- "All miRNA have DHS verified"
    } else {
      fail <- dhs_check$fail_DHS
    }

    if (allmirdhs_byforce == TRUE && !is.null(dhs_check$fail_DHS)) {
      new_mir_tss_df2 <- mir_eponine_score %>%
        filter(mir_name %in% dhs_check$fail_DHS)
      b <- phast_score(new_mir_tss_df2$mir_name, new_mir_tss_df2$chrom,
                       new_mir_tss_df2$stem_loop_p1,
                       new_mir_tss_df2$stem_loop_p2, new_mir_tss_df2$strand,
                       new_mir_tss_df2$tss_p1, new_mir_tss_df2$tss_p2,
                       new_mir_tss_df2$eponine_score)
      a <- bind_rows(a, b)
    }
  }

  list(mir_df = a, fail_DHS = fail)
}



find_candidate_tss_plot <- function(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                                    strand, tss_p1, tss_p2, eponine_score,
                                    ignore_DHS_check = TRUE,
                                    DHS, allmirdhs_byforce = TRUE) {

  mir_eponine_score <- data_frame(mir_name = mir_name, chrom = chrom,
                                  stem_loop_p1 = stem_loop_p1,
                                  stem_loop_p2 = stem_loop_p2,
                                  strand = strand, tss_p1 = tss_p1,
                                  tss_p2 = tss_p2,
                                  eponine_score = eponine_score)


  if (ignore_DHS_check == TRUE) {
    a <- phast_score_plot(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                          strand, tss_p1, tss_p2, eponine_score)

  } else {

    dhs_check <- check_DHS_df(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                              strand, tss_p1, tss_p2, eponine_score, DHS)

    if (is.null(dhs_check$fail_DHS)) {
      a <- phast_score_plot(dhs_check$success$mir_name, dhs_check$success$chrom,
                            dhs_check$success$stem_loop_p1,
                            dhs_check$success$stem_loop_p2,
                            dhs_check$success$strand,
                            dhs_check$success$tss_p1, dhs_check$success$tss_p2,
                            dhs_check$success$eponine_score)
    } else if (allmirdhs_byforce == TRUE) {
      a <- phast_score_plot(mir_name, chrom, stem_loop_p1, stem_loop_p2,
                            strand, tss_p1, tss_p2, eponine_score)
    } else{
      stop("The miRNA you selected failed to pass the DHS check.
           If you want to have the plot anyway,
           you should set the parameter 'allmirdhs_byforce = TRUE'")
    }
  }
  a
}
