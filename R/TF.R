#' Search for putative TFs that regulate miRNA
#'
#' With the putative TSSs of miRNA by \code{\link{find_tss}}, we can further
#' predict transcription factor bind sites by integrating TF motif data from
#' JASPAR2018 and predict transcriptional regulation relationship between TF and
#' miRNA. NOTICE this process may take relatively long time.
#'
#'
#' \code{tss_df} is a data.frame containing as least four columns. These four
#' columns are \code{mir_name}, \code{chrom}, \code{strand} and
#' \code{predicted_tss}.
#'
#' \code{tf_n} means that we will find TFs in upstream tf_n bps from predicted
#' TSSs, which default is 1000.
#'
#' @param tss_df A data.frame. Putative TSSs of miRNA predicted by
#'   \code{\link{find_tss}}.
#' @param tf_n A number. The length of promoter region of miRNA for prediction
#'   of TFBS.
#' @param min.score A single absolute value between 0 and 1. Threshold for
#'   scoring transcription factor binding sites. See also:
#'   \code{\link[TFBSTools]{searchSeq}}.
#'
#' @return A data.frame. The transcriptional regulation relationship between TF
#'   and miRNA.
#'
#' @examples
#'
#' tss_df <- data.frame(
#'                 mir_name = c("hsa-mir-5697", "hsa-mir-192"),
#'                 chrom = c("chr1", "chr11"),
#'                 strand = c("+", "-"),
#'                 predicted_tss = c(10003486, 64683600),
#'                 stringsAsFactors = FALSE)
#' tf <- mir_tf(tss_df, tf_n = 1000)
#'
#' @importFrom TFBSTools getMatrixSet
#' @importFrom TFBSTools toPWM
#' @importFrom Biostrings readDNAStringSet
#' @importFrom TFBSTools searchSeq
#' @importFrom TFBSTools writeGFF3
#' @importFrom stringr str_remove
#' @importFrom tibble as_tibble
#' @importFrom JASPAR2018 JASPAR2018
#'
#' @export


mir_tf <- function(tss_df, tf_n = 1000, min.score = 0.8) {
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
  pfms <- getMatrixSet(JASPAR2018, opts)
  pwms <- toPWM(pfms)
  fasta_read <- readDNAStringSet(tmp_path)

  file.remove(tmp_path)
  sitesetList <- searchSeq(pwms, fasta_read, strand="+", min.score)
  a_GFF <- writeGFF3(sitesetList)
  result <- a_GFF %>%
    separate(attributes, into = c("TF", "TF_class", "TF_sequence"), sep = ";") %>%
    select(seqname, TF, TF_class)

  result$TF <- str_remove(result$TF, "^TF=")
  result$TF_class <- str_remove(result$TF_class, "^class=")

  s_mir_tf <- function(x) {
    unique(x$TF) %>%
      paste(collapse = ", ")
  }

  s_mir_list <- split(result, result$seqname)
  tf <- lapply(s_mir_list, s_mir_tf) %>% unlist()

  data_frame(mir_name = names(s_mir_list), tf = tf)
}
