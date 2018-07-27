#' Predict TSSs of miRNA using graphical web interface.
#'
#' A graphical web interface is provided for users to achieve the functions of
#' \code{\link{find_tss}} and \code{\link{plot_primiRNA}} to intuitively and
#' conveniently predict putative TSSs of miRNA.
#'
#' Users can refer documents of the two functions mentioned ABOVE for details.

#' @examples
#'
#' \dontrun{
#' run_primirTSSapp()
#'}
#'
#' @export
run_primirTSSapp <- function() {
  app_path <- system.file("extdata", "primirTSSapp", package = "primirTSS")
  shiny::runApp(app_path)
}
