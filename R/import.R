#' Import bedgraph files
#'
#' @param filenames Paths to the bedGraph files.
#' @param filter_negative_coverages Convert negative coverage values to 0? (Default: TRUE)
#' @param keep_standard_chromosomes Remove alternative chromosomes? (Default: TRUE)
#'
#' @return A list of GRanges (one element per file).
#'
#' @examples
#' \dontrun{
#'     filenames <- c("file1.bg", "file2.bg")
#'     cov <- import_bedgraphs(filenames)
#' }
#'
#' @import rtracklayer
#' @import purrr
#'
#' @export
import_bedgraphs <- function(filenames,
                             filter_negative_coverages = TRUE,
                             keep_standard_chromosomes = TRUE) {
    stopifnot(all(map_lgl(filenames, file.exists)))
    stopifnot(is.logical(keep_standard_chromosomes))
    stopifnot(is.logical(filter_negative_coverages))

    gr <- map(filenames, import, format = "bedGraph")
    if (filter_negative_coverages) {
        gr <- map(gr, ~ { .x$score[.x$score < 0] <- 0; .x; })
    }
    if (keep_standard_chromosomes) {
        gr <- map(gr, GenomeInfoDb::keepStandardChromosomes,
                  pruning.mode = "coarse")
    }
    gr
}
