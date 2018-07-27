tss_single_judge2 <- function(stem_loop_p1, stem_loop_p2, chrom, strand, tss_p1, tss_p2, gene_loci) {

  gene_df <- gene_loci[gene_loci$chrom == chrom & gene_loci$strand == strand, ]

  intra_index <- stem_loop_p1 >= gene_df$gene_p1 & stem_loop_p2 <= gene_df$gene_p2
  if (any(intra_index)) {
    gene_intra <- gene_df[intra_index, ]
    closest_tss <- lapply(gene_intra$gene_start, function(x) {min(abs(x - tss_p1), abs(x - tss_p2))}) %>%
      unlist()
    in_tss <- (gene_intra$gene_start - tss_p1) * (gene_intra$gene_start - tss_p2) < 0
    colsest_index <- closest_tss <= 150 | in_tss
    if (any(colsest_index)) {
      result <- paste("intra__host_TSS",
                      gene_intra$gene_id[colsest_index],
                      gene_intra$gene_start[colsest_index],
                      # closest_tss[colsest_index],
                      sep = "__") %>%
        paste0(collapse = ":")
    } else {
      tss_site <- ifelse(strand == "+", tss_p2, tss_p1)
      both_site_index <- (tss_site - gene_intra$gene_p1) * (tss_site - gene_intra$gene_p2) < 0
      if (any(both_site_index)) {
        gene_intra_2 <- gene_intra[both_site_index, ]
        distance <- abs(tss_site - gene_intra_2$gene_start)
        gene_intra_2_index <- which.min(distance)
        result <- paste("intra__intra_TSS",
                        gene_intra_2$gene_id[gene_intra_2_index],
                        tss_site,
                        sep = "__")
      } else {
        distance <- abs(tss_site - gene_intra$gene_start)
        gene_intra_2_index <- which.min(distance)
        result <- paste("intra__intra_TSS2",
                        gene_intra$gene_id[gene_intra_2_index],
                        tss_site,
                        sep = "__")
      }
    }
  } else {
    closest_tss <- lapply(gene_df$gene_start, function(x) {min(abs(x - tss_p1), abs(x - tss_p2))}) %>%
      unlist()
    in_tss <- (gene_df$gene_start - tss_p1) * (gene_df$gene_start - tss_p2) < 0
    colsest_index <- closest_tss <= 150 | in_tss

    if (any(colsest_index)) {
      tss_site <- ifelse(strand == "+", tss_p2, tss_p1)
      mir_end <- ifelse(strand == "+", stem_loop_p2, stem_loop_p1)

      distance <- abs(tss_site - gene_df$gene_start)
      gene_inter_index <- which.min(distance)
      result <- paste("inter__overlap_inter_TSS2",
                      gene_df$gene_id[gene_inter_index],
                      tss_site,
                      sep = "__")
    } else {
      tss_site <- ifelse(strand == "+", tss_p2, tss_p1)
      mir_end <- ifelse(strand == "+", stem_loop_p2, stem_loop_p1)
      gene_tss_trans_index <- (tss_site - gene_df$gene_start) * (mir_end - gene_df$gene_start) < 0

      if (any(gene_tss_trans_index)) {
        distance <- abs(tss_site - gene_df$gene_start)
        gene_inter_index <- which.min(distance)
        result <- paste("inter__overlap_inter_TSS2",
                        gene_df$gene_id[gene_inter_index],
                        tss_site,
                        sep = "__")
      } else {
        gene_end <- ifelse(gene_df$strand == "+", gene_df$gene_p2, gene_df$gene_p1)
        gene_end_trans_index <- (tss_site - gene_end) * (mir_end - gene_end) < 0
        if (any(gene_end_trans_index)) {
          select_gene <- gene_df[gene_end_trans_index, ]
          distance <- abs(tss_site - select_gene$gene_start)
          gene_inter_index <- which.min(distance)
          result <- paste("inter__overlap_inter_TSS",
                          select_gene$gene_id[gene_inter_index],
                          tss_site,
                          sep = "__")
        } else {
          distance <- abs(tss_site - gene_df$gene_start)
          gene_inter_index <- which.min(distance)
          result <- paste("inter__inter_TSS",
                          gene_df$gene_id[gene_inter_index],
                          tss_site,
                          sep = "__")
        }
      }
    }
  }
  result
}


tss_filter2 <- function(candidate_tss, expressed_gene = "all") {


  if (length(expressed_gene) > 1 || expressed_gene != "all") {
    gene_loci <- gene_loci %>%
      filter(gene_id %in% expressed_gene)
  }

  tss_info <- candidate_tss %>%
    select(stem_loop_p1, stem_loop_p2, chrom, strand, tss_p1, tss_p2) %>%
    pmap(tss_single_judge2, gene_loci) %>%
    unlist()

  s_tss <- candidate_tss %>%
    bind_cols(tss_info = tss_info) %>%
    select(-tss_p1, -tss_p2) %>%
    filter(tss_info != "*") %>%
    extand_df("tss_info", ":") %>%
    separate(new_info,
             into = c("mir_context", "tss_type", "gene", "predicted_tss"),
             sep = "__") %>%
    mutate(predicted_tss = as.double(predicted_tss))

  a1 <- s_tss$stem_loop_p1 - s_tss$predicted_tss
  a2 <- s_tss$predicted_tss - s_tss$stem_loop_p2
  pri_tss_distance <- ifelse(s_tss$strand == "+", a1, a2)
  bind_cols(s_tss, pri_tss_distance = pri_tss_distance) %>%
    group_by(mir_name) %>%
    filter(pri_tss_distance == min(pri_tss_distance)) %>%
    ungroup()
}















