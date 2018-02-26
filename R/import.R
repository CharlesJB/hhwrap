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
#' @param filename Path to the peak file.
#' @param genome Add seqinfo from genome (see ?fetchExtendedChromInfoFromUCSC
#'        for the list of available genomes). Value must be NULL or a character
#'        string (i.e.: "hg38").
#' @param keep_standard_chromosomes Remove alternative chromosomes? (Default: TRUE)
#'
#' @return A GRanges corresponding the regions in the peak file.
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
import_peaks <- function(filename,
                         genome = NULL,
                         keep_standard_chromosomes = TRUE) {

    stopifnot(length(filename) == 1)
    stopifnot(file.exists(filename))
    stopifnot(str_detect(filename, "\\.(narrow|broad)Peak"))

    if (!is.null(genome)) {
        data(si)
        stopifnot(is.character(genome))
        stopifnot(length(genome) == 1)
        stopifnot(any(names(si) == genome))
        genome <- si[[genome]]
        if (keep_standard_chromosomes) {
            genome <- GenomeInfoDb::keepStandardChromosomes(genome)
        }
    }

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
    extraCols <- get_extraCols(str_extract(filename, "narrowPeak$|broadPeak$"))

    if (!is.null(genome)) {
        peak <- import(filename, extraCols = extraCols, format = "BED",
                       genome = genome)
    } else {
        peak <- import(filename, extraCols = extraCols, format = "BED")
    }
    if (keep_standard_chromosomes) {
        peak <- GenomeInfoDb::keepStandardChromosomes(peak)
    }
    peak
}
