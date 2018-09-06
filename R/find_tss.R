#' Predict TSSs of miRNA
#'
#' Search for putative TSSs of miRNA, together with integrating available data such
#' as H3K4me3 data, Pol II data, miRNA expression data, and protein-coding gene
#' data, as well as provide the transcriptional regulation relationship between TF
#' and miRNA.
#'
#' @param bed_merged Peaks from ChIP-seq data to be provided for analysis can be
#'   H3K4me3 peaks, Pol II peaks or both. Notice that peaks are supposed to be
#'   merged(see also \code{\link{peak_merge}}) before \code{find_TSS} if using
#'   only one kind of peak data, while peaks should be firstly merged and then
#'   join together(see also \code{\link{peak_join}}) if both H3K4me3 data and
#'   Pol II are input.
#' @param expressed_mir This parameter allows users to specify certain miRNAs,
#'   the TSSs of which they want to search for by providing a list of
#'   miRNAs(e.g., expressed miRNAs in a certain cell-line). If
#'   \code{expressed_mir} is not specified, the default value of the parameter
#'   is "\code{all}" and the function will acquiescently employ all the miRNAs
#'   currently listed on "\code{miRbase}" database.
#' @param flanking_num A parameter in Eponine model to detect TSSs. It is
#'   concluded that a peak signal with flanking regions of C-G enrichment are
#'   important to mark TSSs. The default value is 1000.
#' @param threshold The threshold for candidate TSSs scored with Eponine method. The
#'   default value is 0.7.
#' @param ignore_DHS_check The process of DHS_check further assists to filter
#'   putative TSSs. When there is a DHS peak that locates within 1 kb upstream
#'   of a putative TSS, this predicted TSS will be retained for its character is
#'   consistent with that of an authentic TSS. Or the TSSs with no DHSs locating
#'   within 1 kb upstream of them would be discarded.
#' @param DHS ChIP-seq data of DNase I hypersensitive sites(DHSs).
#' @param allmirdhs_byforce When we use DHS data to check the validity of TSSs,
#'   there is a possibility where no DHSs locates within 1 kb upstream of all
#'   putative TSSs and all these putative TSSs might be filtered out by our
#'   method resulting no outputs. While "\code{allmirdhs_byforce} = TRUE", it
#'   ensures to output at least 1 most possible TSS even if the nearest DHS
#'   signal locates more than 1 kb upstream of this TSS.
#' @param expressed_gene Users can specify genes expressed in certain
#'   cell-lines that are analyzed. Or the default value is "\code{all}", which
#'   means all the expressed genes annotated on Ensemble will be employed.
#' @param allmirgene_byforce While integrating expressed_gene data to improve
#'   prediction, there might be a circumstance where all the putative TSS are
#'   discarded. To prevent this condition, users are allowed to use
#'   "\code{allmirgene_byforce = TRUE}" to ensure at least 1 putative TSS for
#'   each miRNA will be output.
#' @param seek_tf With the result of predicted TSSs, seek_tf provides users with
#'   an option to predict related TFs for miRNA. The data of transcription
#'   factors refer to \code{JASPAR2018} database.
#' @param tf_n TFBS locates on the upstream of the TSS of a certain TF, which is
#'   considered as the promoter region. \code{tf_n} set the length of promoter
#'   region for predicting transcription regulation between miRNAs and TFs.
#' @param min.score The threshold for scoring transcription factor binding sites. A
#'   single absolute value between 0 and 1.
#'
#'
#'
#' @return The first part of the result returns details of predicted TSSs,
#'   composed of seven columns: \code{mir_name, chrom, stem_loop_p1,
#'   stem_loop_p2, strand mir_context, tss_type gene and predicted_tss}:
#'
#'   \code{mir_name}: Name of miRNA.
#'
#'   \code{chrom}: Chromosome.
#'
#'   \code{stem_loop_p1}: The start site of a stem-loop.
#'
#'   \code{stem_loop_p2}: The end site of a stem-loop.
#'
#'   \code{strand}: Polynucleotide strands. (\code{+/-})
#'
#'   \code{mir_context}: The relative positon relationship between stem-loop and
#'   protein-coding gene. (\code{intra/inter})
#'
#'   \code{tss_type}: Four types of predicted TSSs. See the section below
#'   \code{TSS types} for details.
#'   (\code{host_TSS/intra_TSS/overlap_inter_TSS/inter_TSS})
#'
#'   \code{gene}: Ensembl gene ID
#'
#'   \code{predicted_tss}: Predicted transcription start sites(TSSs).
#'
#'   \code{pri_tss_distance}: The distance between a predicted TSS and the start
#'   site of the stem-loop.
#'
#'
#' @section TSS types: TSSs are catalogued into 4 types as below.
#'
#'   \code{host_TSS} The TSSs of miRNA that are close to the TSS of
#'   protein-coding gene implying they may share the same TSS, on the condition
#'   where \code{mir_context} is "\code{intra}". (See above:
#'   \code{Value-mir_context})
#'
#'   \code{intra_TSS} The TSSs of miRNA that are NOT close to the TSS of
#'   the protein-coding gene, on the condition where \code{mir_context} is
#'   "\code{intra}".
#'
#'   \code{overlap_inter_TSS} The TSSs of miRNA are catalogued as
#'   "\code{overlap_inter_TSS}" when the pri-miRNA gene overlaps with Ensembl
#'   gene, on the condition where "\code{mir_context}" is "\code{inter}".
#'
#'   \code{inter_inter_TSS} The TSSs of miRNA are catalogued as
#'   "\code{inter_inter_TSS}" when the miRNA gene does NOT overlap with Ensembl
#'   gene, on the condition where "\code{mir_context}" is "\code{inter}".
#'
#'   (See Xu HUA et al 2016 for more details)
#'
#'
#'
#' @section Log:
#'
#'   The second part of the result returns logs during the process of
#'   prediction: \code{find_nearest_peak_log} If no peaks locate in the
#'   upstream of a stem-loop to help determine putative TSSs of miRNA, we will
#'   fail to find the nearest peak and this miRNA will be logged in
#'   \code{find_nearest_peak_log}.
#'
#'   \code{eponine_score_log} For a certain miRNA, if none of the candidate
#'   TSSs scored with \code{Eponine method} meet the threshold we set, we will
#'   fail to get an eponine score and this miRNA will be logged in
#'   \code{eponine_score_log}.
#'
#'   \code{DHS_check_log} For a certain miRNA, if no DHS signals locate
#'   within 1 kb upstream of each putative TSSs, these putative TSSs will be
#'   filtered out and this miRNA will be logged in \code{DHS_check_log}.
#'
#'   \code{gene_filter_log} For a certain miRNA, when integrating
#'   expressed_gene data to improve prediction, if no putative TSSs are
#'   confirmed after considering the relative position relationship among TSSs,
#'   stem-loops and expressed genes, this miRNA will be filtered out and logged
#'   in \code{gene_filter_log}.
#'
#' @section Reference: Xu Hua, Luxiao Chen, Jin Wang*, Jie Li* and Edgar
#'   Wingender*, Identifying cell-specific microRNA transcriptional start sites.
#'   Bioinformatics 2016, 32(16), 2403-10.
#'
#'
#' @examples
#'
#' bed_merged <- data.frame(
#'                 chrom = c("chr1", "chr1", "chr1", "chr1", "chr2"),
#'                 start = c(9910686, 9942202, 9996940, 10032962, 9830615),
#'                 end = c(9911113, 9944469, 9998065, 10035458, 9917994),
#'                 stringsAsFactors = FALSE)
#' bed_merged <- as(bed_merged, "GRanges")
#'
#' ownmiRNA <- find_tss(bed_merged, expressed_mir = "hsa-mir-5697",
#'                      ignore_DHS_check = TRUE,
#'                      expressed_gene = "all",
#'                      allmirgene_byforce = TRUE)
#'
#'
#' @export
find_tss <- function(bed_merged, expressed_mir = "all",
                     flanking_num = 1000, threshold = 0.7,
                     ignore_DHS_check = TRUE,
                     DHS, allmirdhs_byforce = TRUE,
                     expressed_gene = "all", allmirgene_byforce = TRUE,
                     seek_tf = FALSE, tf_n = 1000, min.score = 0.8) {

  h3_mir_flank <- find_nearest_peak(bed_merged, expressed_mir)

  mir_eponine_score <- eponine_score(h3_mir_flank$success$mir_name,
                                     h3_mir_flank$success$chrom,
                                     h3_mir_flank$success$stem_loop_p1,
                                     h3_mir_flank$success$stem_loop_p2,
                                     h3_mir_flank$success$strand,
                                     h3_mir_flank$success$peak_p1,
                                     h3_mir_flank$success$peak_p2,
                                     flanking_num, threshold)

  candidate_tss <- find_candidate_tss(mir_eponine_score$success$mir_name,
                                      mir_eponine_score$success$chrom,
                                      mir_eponine_score$success$stem_loop_p1,
                                      mir_eponine_score$success$stem_loop_p2,
                                      mir_eponine_score$success$strand,
                                      mir_eponine_score$success$tss_p1,
                                      mir_eponine_score$success$tss_p2,
                                      mir_eponine_score$success$eponine_score,
                                      ignore_DHS_check,
                                      DHS, allmirdhs_byforce)

  gene_fliter_tss <- tss_filter(candidate_tss$mir_df$mir_name,
                                candidate_tss$mir_df$chrom,
                                candidate_tss$mir_df$stem_loop_p1,
                                candidate_tss$mir_df$stem_loop_p2,
                                candidate_tss$mir_df$strand,
                                candidate_tss$mir_df$tss_p1,
                                candidate_tss$mir_df$tss_p2,
                                expressed_gene, allmirgene_byforce)

  if (any(h3_mir_flank$fail_nearpeak !=
          "No one miRNA fail to find its nearest peaks")) {
    if (length(h3_mir_flank$fail_nearpeak) == 1) {
      log_1 <- paste(h3_mir_flank$fail_nearpeak,
                     "fail to find its nearest peaks", " ")
    } else {
      log_1 <- paste(h3_mir_flank$fail_nearpeak, collapse = ", ")
      log_1 <- paste(log_1, "fail to find their nearest peaks", " ")
    }
  } else {
    log_1 <- h3_mir_flank$fail_nearpeak
  }

  if (any(mir_eponine_score$fail_eponine != "All miRNAs have eponine socres")) {
    if (length(mir_eponine_score$fail_eponine) == 1) {
      log_2 <- paste(mir_eponine_score$fail_eponine,
                     "fail to have eponine socres", " ")
    } else {
      log_2 <- paste(mir_eponine_score$fail_eponine, collapse = ", ")
      log_2 <- paste(log_2, "fail to have eponine socres", " ")
    }
  } else {
    log_2 <- mir_eponine_score$fail_eponine
  }

  if (any(candidate_tss$fail_DHS != "Do not have a DHS check") &&
      any(candidate_tss$fail_DHS != "All miRNA have DHS verified")) {
    if (length(candidate_tss$fail_DHS) == 1) {
      log_3 <- paste(candidate_tss$fail_DHS, "fail to have a DHS check", " ")
    } else {
      log_3 <- paste(candidate_tss$fail_DHS, collapse = ", ")
      log_3 <- paste(log_3, "fail to have a DHS check", " ")
    }
  } else {
    log_3 <- candidate_tss$fail_DHS
  }

  if (any(gene_fliter_tss$fail_genefilter != "All miRNAs have tss")) {
    if (length(gene_fliter_tss$fail_genefilter) == 1) {
      log_4 <- paste(gene_fliter_tss$fail_genefilter, "fail to have a tss", " ")
    } else {
      log_4 <- paste(gene_fliter_tss$fail_genefilter, collapse = ", ")
      log_4 <- paste(log_4, "fail to have a tss", " ")
    }
  } else {
    log_4 <- gene_fliter_tss$fail_genefilter
  }

  if (allmirgene_byforce == TRUE) {
    log_4 <- "All miRNAs have tss"
  }

  if (seek_tf == TRUE) {
    a <- mir_tf(gene_fliter_tss$tss_df$mir_name, gene_fliter_tss$tss_df$chrom,
                gene_fliter_tss$tss_df$strand,
                gene_fliter_tss$tss_df$predicted_tss,
                tf_n, min.score)

    result <- dplyr::left_join(gene_fliter_tss$tss_df, a, by = "mir_name")
  } else {
    result <- gene_fliter_tss$tss_df
  }

  list(tss_df = result,
       find_nearest_peak_log = log_1,
       eponine_score_log = log_2,
       DHS_check_log = log_3,
       gene_filter_log = log_4)
}






