#' Split a GRanges into N bins
#'
#' Originally posted by Martin Morgan:
#' https://stat.ethz.ch/pipermail/bioconductor/2012-September/047923.html
#'
#' Wparam gr A GRanges with only one seqnames value.
#' @param n Number of bins to produce.
#'
#' @return
#'   A GRanges object splitted into N bins
#'
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges width
#' @importFrom IRanges IRanges
#'
#' @examples
#'   gr <- GRanges("chr1", IRanges(c(100, 300), c(200, 500))
#'   gr <- intoNbins(gr)
intoNbins <- function(gr, n = 10) {
    if (any(width(gr) < n)) {
        stop("all 'width(gr)' must be >= 'n'")
    }
    d <- width(gr) / n
    dd <- cumsum(rep(d, each=n))
    mask <- logical(n); mask[1] <- TRUE
    dd <- dd - rep(dd[mask], each=n)

    starts <- round(rep(start(gr), each=n) + dd)
    ends <- c(starts[-1], 0) - 1L
    ends[rev(mask)] <- end(gr)

    gr <- gr[rep(seq_along(gr), each=n)]
    GenomicRanges::ranges(gr) <- IRanges(starts, ends)
    gr
}

#' Convert coverages to matrix
#'
#' @param coverage A single coverage obtained with the import_bedgraphs
#'                 function.
#' @param peaks A GRanges corresponding to the subset of regions to subset
#'              the coverages.
#' @param ncol The number of columns in the result matrix. (Default = 100).
#'
#' @return A matrix with ncol columns and length(peaks) rows.
#'
#' @examples
#' filenames <- get_demo_bedGraph_files()
#' cov <- import_bedgraphs(filenames)
#' 
#' peaks <- import_peaks(get_demo_peaks_file())
#' m <- coverage_2_matrix(cov[[1]], peaks)
#'
#' @import purrr
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges Views
#' @importFrom IRanges viewMeans
#' @import GenomeInfoDb
#'
#' @export
coverage_2_matrix <- function(coverage, peaks, ncol = 100) {
    stopifnot(is(peaks, "GRanges"))
    stopifnot(all(GenomicRanges::width(peaks) >= ncol))

    i <- S4Vectors::queryHits(GenomicRanges::findOverlaps(coverage, peaks))
    coverage <- coverage[i]
    coverage <- GenomicRanges::coverage(coverage, weight = coverage$score)

    peaks <- split(peaks, as.character(seqnames(peaks)))

    extract_scores <- function(n) {
        gr <- intoNbins(peaks[[n]], n = ncol)
        cov <- coverage[[n]]
        view <- Views(cov, start(gr), end(gr))
        matrix(viewMeans(view), ncol = ncol, byrow = TRUE)
    }
    m <- map(names(peaks), extract_scores)
    do.call("rbind", m)
}
