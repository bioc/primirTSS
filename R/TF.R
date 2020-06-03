#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_remove
#' @importFrom tibble as_tibble
#' @importFrom JASPAR2018 JASPAR2018
#'


mir_tf <- function(mir_name, chrom, strand, predicted_tss,
                   tf_n = 1000, min.score = 0.8) {

  tss_df <- tibble::tibble(mir_name = mir_name,
                       chrom = chrom,
                       strand = strand,
                       predicted_tss = predicted_tss)

  a <- tss_df %>%
    select(mir_name, chrom, strand, predicted_tss)

  a1 <- paste(a$predicted_tss - tf_n, a$predicted_tss, sep = "--")
  a2 <- paste(a$predicted_tss, a$predicted_tss + tf_n, sep = "--")

  region <- ifelse(a$strand == "+", a1, a2)
  tcga <- bind_cols(a, region = region) %>%
    separate(region, into = c("start", "end")) %>%
    select(-predicted_tss) %>%
    mutate(start = as.double(start),
           end = as.double(end)) %>%
    select(chrom, start, end, strand) %>%
    pmap(cor_tdga) %>%
    unlist()

  line1 <- a %>%
    mutate(arrow = ">") %>%
    unite(line1, arrow, mir_name, sep = "")

  a_tmp <- paste(line1$line1, tcga, sep = "\n")
  tmp_path <- file.path(tempdir(), paste0("tf", runif(1, 0, 100), ".fa"))
  writeLines(a_tmp, tmp_path)

  opts <- list()
  opts[["species"]] <- "9606"
  pfms <- TFBSTools::getMatrixSet(JASPAR2018, opts)
  pwms <- TFBSTools::toPWM(pfms)
  fasta_read <- readDNAStringSet(tmp_path)

  file.remove(tmp_path)
  sitesetList <- TFBSTools::searchSeq(pwms, fasta_read, strand="+", min.score)
  a_GFF <- TFBSTools::writeGFF3(sitesetList)
  result <- a_GFF %>%
    separate(attributes,
             into = c("TF", "TF_class", "TF_sequence"),
             sep = ";") %>%
    select(seqname, TF, TF_class)

  result$TF <- str_remove(result$TF, "^TF=")
  result$TF_class <- str_remove(result$TF_class, "^class=")

  s_mir_tf <- function(x) {
    unique(x$TF) %>%
      paste(collapse = ", ")
  }

  s_mir_list <- split(result, result$seqname)
  tf <- lapply(s_mir_list, s_mir_tf) %>% unlist()

  tibble::tibble(mir_name = names(s_mir_list), tf = tf)
}
