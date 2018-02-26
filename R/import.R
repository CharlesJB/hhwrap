#' Import bedgraph files
#'
#' @param filenames Paths to the bedGraph files.
#' @param filter_negative_coverages Convert negative coverage values to 0? (Default: TRUE)
#' @param keep_standard_chromosomes Remove alternative chromosomes? (Default: TRUE)
#'
#' @return A list of GRanges (one element per file).
#'
#' @examples
#' filenames <- get_demo_bedGraph_files()
#' cov <- import_bedgraphs(filenames)
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

#' Import peak files
#'
#' Only works for .narrowPeak or .broadPeak file formats.
#'
#' @param filenames Paths to the peak files.
#'
#' @return A list of GRanges (one element per file).
#'
#' @examples
#' peaks_file <- get_demo_peaks_file()
#' peaks <- import_peaks(peaks_file)
#'
#' @import rtracklayer
#' @import stringr
#' @import tools
#'
#' @export
import_peaks <- function(filenames) {
    stopifnot(all(str_detect(filenames, "\\.(narrow|broad)Peak")))
    stopifnot(all(map_lgl(filenames, file.exists)))

    get_extraCols <- function(x) {
        if (x == "narrowPeak") {
            c(signalValue = "numeric",
              pValue = "numeric",
              qValue = "numeric",
              peak = "integer")
        } else {
            c(signalValue = "numeric",
              pValue = "numeric",
              qValue = "numeric")
        }
    }

    extraCols <- map(str_extract(filenames, "narrowPeak$|broadPeak$"),
                     get_extraCols)

    names(filenames) <- file_path_sans_ext(basename(filenames))
    map2(filenames, extraCols, ~ import(.x, extraCols = .y, format = "BED"))
}
