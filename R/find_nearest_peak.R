

find_nearest_peak <- function(peak, expressed_mir = "all") {
  peak <- as.data.frame(peak)
  peak <- peak %>%
    mutate(chrom = as.character(seqnames),
           start = as.double(start),
           end = as.double(end)) %>%
    select(chrom, start, end)

  peak_list <- split(peak, peak$chrom)

  if (length(expressed_mir) > 1 || expressed_mir != "all") {
    id_mir <- mir_loci$mir_name
    mir_loci <- mir_loci[id_mir %in% expressed_mir, ]
  }

  mir_loci_list <- split(mir_loci, mir_loci$chrom)

  all_chr <- intersect(names(mir_loci_list), names(peak_list))

  split_result <- vector("list", length(all_chr))

  for (i in seq_along(all_chr)) {
    tmp_peak <- peak_list[[all_chr[i]]]
    tmp_mir <- mir_loci_list[[all_chr[i]]]

    f_fun <- function(x) {
      compare <- x - tmp_peak$end
      if(all(compare >= 0)) return("*")
      compare[compare >= 0] <- -Inf
      which.max(compare) %>% as.character()
    }
    aa1 <- lapply(tmp_mir$mir_start, f_fun) %>% unlist()

    z_fun <- function(x) {
      compare <- x - tmp_peak$start
      if(all(compare <= 0)) return("*")
      compare[compare <= 0] <- Inf
      which.min(compare) %>% as.character()
    }
    aa2 <- lapply(tmp_mir$mir_start, z_fun) %>% unlist()

    index <- ifelse(tmp_mir$strand == "-", aa1, aa2)

    if (any(index == "*")) {
      rm_index <- which(index == "*")
      tmp_mir <- tmp_mir[-rm_index, ]
      index <- index[-rm_index]
    }

    index <- as.double(index)

    tmp_df <- bind_cols(tmp_mir, tmp_peak[index, ])
    index1 <-
      ((tmp_df$mir_start - tmp_df$start1) * (tmp_df$mir_start - tmp_df$end1)) <=
      0
    tmp_df[(index1 & tmp_df$strand == "+"), ]$end1 <-
      tmp_df[(index1 & tmp_df$strand == "+"), ]$mir_start - 1
    tmp_df[(index1 & tmp_df$strand == "-"), ]$start1 <-
      tmp_df[(index1 & tmp_df$strand == "-"), ]$mir_start + 1

    split_result[[i]] <- tmp_df %>%
      select(mir_name, chrom, stem_loop_p1 = start, stem_loop_p2 = end, strand,
             peak_p1 = start1, peak_p2 = end1)
  }

  success_tf <- Reduce(rbind, split_result)

  fail_res <- setdiff(unique(mir_loci$mir_name), unique(success_tf$mir_name))
  if (length(fail_res) == 0) {
    fail <- "No one miRNA fail to find its nearest peaks"
  } else {
    fail <- fail_res
  }
  list(success = Reduce(rbind, split_result), fail_nearpeak = fail)
}
