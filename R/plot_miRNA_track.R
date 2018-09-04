#' @importFrom dplyr bind_rows
#' @importFrom Gviz GeneRegionTrack
#' @importFrom Gviz DataTrack
#'

plot_primiRNA_track <- function(expressed_mir, bed_merged,
                          flanking_num = 1000, threshold = 0.7,
                          ignore_DHS_check = TRUE,
                          DHS, allmirdhs_byforce = TRUE,
                          expressed_gene = "all",
                          allmirgene_byforce = TRUE) {

  genome_version <- "hg38"
  mir_peaks <- find_nearest_peak(bed_merged, expressed_mir)
  mir_eponine_score <- eponine_score(mir_peaks$success$mir_name,
                                     mir_peaks$success$chrom,
                                     mir_peaks$success$stem_loop_p1,
                                     mir_peaks$success$stem_loop_p2,
                                     mir_peaks$success$strand,
                                     mir_peaks$success$peak_p1,
                                     mir_peaks$success$peak_p2,
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

  tss_list <- tss_filter(candidate_tss$mir_df$mir_name,
                         candidate_tss$mir_df$chrom,
                         candidate_tss$mir_df$stem_loop_p1,
                         candidate_tss$mir_df$stem_loop_p2,
                         candidate_tss$mir_df$strand,
                         candidate_tss$mir_df$tss_p1,
                         candidate_tss$mir_df$tss_p2,
                         expressed_gene, allmirgene_byforce)

  tss <- tss_list$tss_df

  ep_con <- find_candidate_tss_plot(mir_eponine_score$success$mir_name,
                                    mir_eponine_score$success$chrom,
                                    mir_eponine_score$success$stem_loop_p1,
                                    mir_eponine_score$success$stem_loop_p2,
                                    mir_eponine_score$success$strand,
                                    mir_eponine_score$success$tss_p1,
                                    mir_eponine_score$success$tss_p2,
                                    mir_eponine_score$success$eponine_score,
                                    ignore_DHS_check,
                                    DHS, allmirdhs_byforce)

  tss_tmp <- tss %>%
    mutate(start = predicted_tss,
           end = predicted_tss + 1,
           symbol = "tss") %>%
    select(chrom, start, end, strand, symbol)
  tss_p <- tss %>%
    mutate(symbol = "stem-loop") %>%
    select(chrom, start = stem_loop_p1, end = stem_loop_p2, strand, symbol) %>%
    bind_rows(tss_tmp)

  chr <- unique(tss_p$chrom)
  tsstrack <- GeneRegionTrack(tss_p, genome = genome_version,
                             chromosome = chr, name = "pri-miRNA",
                             showId = TRUE,
                             transcriptAnnotation = "symbol",
                             just.group = "above",
                             showOverplotting = FALSE)

  gene_p <- gene_loci[gene_loci$gene_id == tss$gene, ] %>%
    mutate(symbol = gene_id) %>%
    select(start = gene_p1, end = gene_p2, strand, symbol)
  genetrack <- GeneRegionTrack(gene_p, genome = genome_version,
                               chromosome = chr, name = "Ensemble genes",
                               shape = "arrow",
                               collapseTranscripts = "meta",
                               just.group = "above",
                               transcriptAnnotation = "symbol")

  eponine_p <- ep_con %>%
    mutate(start = tss_p1,
           end = tss_p2) %>%
    select(start, end, strand, data = eponine_score)

  e_track = DataTrack(range = eponine_p, genome = genome_version,
                      chromosome = chr,
                      name = "eponine score",
                      type = "h",
                      col.line = "skyblue",
                      ylim = c(0, 1.2))

  con_p <- ep_con %>%
    mutate(start = tss_p1,
           end = tss_p2) %>%
    select(start, end, strand, data = phast_score)

  con_track = DataTrack(range = con_p, genome = genome_version,
                      chromosome = chr,
                      name = "conservation score",
                      type = "h",
                      col.line = "skyblue",
                      ylim = c(0, 1.2))

  min_loci <- min(gene_p$start, gene_p$end, tss_p$start, tss_p$end)
  max_loci <- max(gene_p$start, gene_p$end, tss_p$start, tss_p$end)

  list(tsstrack = tsstrack,
       genetrack = genetrack,
       e_track = e_track,
       con_track = con_track,
       min_loci = min_loci,
       max_loci= max_loci,
       chr = chr)
}


#' Plot the result of prediction for miRNA
#'
#' For each miRNA, plot the position of TSS, pri-miRNA, related Ensemble gene,
#' eponine socre and conservation score according to the result of prediction
#' using primirTSS.
#'
#' NOTICE that this function is used for visualizing the predicted result of
#' ONLY ONE specific miRNA every single time.
#'
#' @param expressed_mir This parameter allows users to specify certain miRNAs,
#'   the TSSs of which they want to search for by providing a list of
#'   miRNAs(e.g. expressed miRNAs in a certain cell-line). If
#'   \code{expressed_mir} is not specified, the default value of the parameter
#'   is "\code{all}" and the function will acquiescently employ all the miRNAs
#'   currently listed on "\code{miRbase}" database.
#' @param bed_merged Peaks from ChIP-seq data to be provided for analysis can be
#'   H3K4me3 peaks, Pol II peaks or both. Notice that peaks are supposed to be
#'   merged(see also \code{\link{peak_merge}}) before \code{find_TSS} if using
#'   only one kind of peak data, while peaks shoud be firstly merged and then
#'   join together(see also \code{\link{peak_join}}) if both H3K4me3 data and
#'   Pol II are input.
#' @param flanking_num A parameter in Eponine model to detect TSSs. It is
#'   concluded that a peak signal with flanking regions of C-G enrichment are
#'   important to mark TSSs. The default value is 1000.
#' @param threshold Threshold for candidate TSSs scored with Eponine method. The
#'   default value is 0.7.
#' @param ignore_DHS_check The process of DHS_check further assist to filter
#'   putative TSSs. When there are a DHS peak that locates within 1 kb upstream
#'   of a putative TSS, this predicted TSS will be retain for it character is
#'   consistent with that of an authentic TSS. Or the TSSs with no DHSs locating
#'   within 1 kb upstream of them would be discard.
#' @param DHS ChIP-seq data of DNase I hypersensitive sites(DHSs).
#' @param allmirdhs_byforce When we use DHS data to check the validity of TSSs,
#'   there is possibility where no DHSs locates within 1 kb upstream of all
#'   putative TSSs and all these putative TSSs might be filtered out by our
#'   method resulting no outputs. While "\code{allmirdhs_byforce} = TRUE", it
#'   ensures to output at least 1 most possible TSS even if the nearest DHS
#'   signal locates more than 1 kb upsteam of this TSS.
#' @param expressed_gene Users can speicify genes expressed in certain
#'   cell-lines that is analyzed. Or the default value is "\code{all}", which
#'   means all the expressed genes annotated on Ensemble will be employed.
#' @param allmirgene_byforce While integrating expressed_gene data to improve
#'   prediction, there might be a circumstance where all the putative TSS are
#'   discarded. To prevent this condition, users are allowed to use
#'   "\code{allmirgene_byforce = TRUE}" to ensure at least 1 putative TSS for
#'   each miRNA will be output.
#'
#' @return There will be six tracks plotted as return:
#'
#' \code{Chrom}: Position of miRNA on the chromosome.
#'
#' \code{hg38}: Reference genome coordinate in hg38.
#'
#' \code{pri-miRNA}: Position of pri-miRNA.
#'
#' \code{Ensemble genes}: Position of related protein-coding gene.
#'
#' \code{eponine score}: Score of best putative TSS conducted by eponine method.
#'
#' \code{conservation score}: Conservation score should be integrated with
#' eponine score to find out putative TSSs.
#'
#'
#'
#'
#' @examples
#' expressed_mir <- "hsa-mir-5697"
#' bed_merged <- data.frame(
#'   chrom = c("chr1", "chr1", "chr1", "chr1", "chr2"),
#'   start = c(9180799, 9201483, 9234339, 9942202, 9830615),
#'   end = c(9183889, 9202580, 9235853, 9944469, 9917994),
#'   stringsAsFactors = FALSE
#' )
#' bed_merged <- as(bed_merged, "GRanges")
#' \donttest{
#' plot_primiRNA(expressed_mir, bed_merged)
#' }
#'
#' @importFrom Gviz IdeogramTrack
#' @importFrom Gviz GenomeAxisTrack
#' @importFrom Gviz plotTracks
#' @export

plot_primiRNA <- function(expressed_mir, bed_merged,
                          flanking_num = 1000, threshold = 0.7,
                          ignore_DHS_check = TRUE,
                          DHS, allmirdhs_byforce = TRUE,
                          expressed_gene = "all",
                          allmirgene_byforce = TRUE) {

  track <- plot_primiRNA_track(expressed_mir, bed_merged,
                               flanking_num, threshold,
                               ignore_DHS_check,
                               DHS, allmirdhs_byforce,
                               expressed_gene,
                               allmirgene_byforce)

  genome_version <- "hg38"
  itrack <- IdeogramTrack(genome = genome_version, showBandId = TRUE, name = "")
  axistrack <- GenomeAxisTrack(add53=TRUE, exponent=0,
                               littleTicks=TRUE, name=genome_version)

  tsstrack <- track$tsstrack
  genetrack <- track$genetrack
  e_track <- track$e_track
  con_track <- track$con_track
  chr <- track$chr
  START <- track$min_loci
  END <- track$max_loci

  plotTracks(
    list(itrack,
         axistrack,
         tsstrack,
         genetrack,
         e_track,
         con_track),
    chromosome = chr,
    from = START,
    to = END,
    extend.left = 0.1,
    extend.right = 0.1,
    fontcolor.group = "black",
    fontcolor = "black",
    background.title = "darkblue",
    col = NULL,
    showTitle = TRUE,
    frame = TRUE,
    fontsize = 11,
    fontsize.group = 11,
    cex = 0.8,
    cex.group = 0.8,
    cex.title = 1,
    sizes = c(1.5, 1.5, 1.5, 1.5, 3, 3)
  )
}

