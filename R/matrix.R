#' Split a GRanges into N bins
#'
#' Originally posted by Martin Morgan:
#' https://stat.ethz.ch/pipermail/bioconductor/2012-September/047923.html
#'
#' @param gr A GRanges with only one seqnames value.
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
#' @param gr A GRanges corresponding to the regions to subset the coverages.
#'           All the regions must have the same width (see
#'           ?GenomicRanges::resize to resize regions).
#' @param ncol The number of columns in the result matrix. (Default = 100).
#'
#' @return A matrix with ncol columns and length(gr) rows.
#'
#' @examples
#' filenames <- get_demo_bedGraph_files()
#' cov <- import_bedgraphs(filenames)
#' 
#' gr <- import_peaks(get_demo_peaks_file())
#' m <- coverage_2_matrix(cov[[1]], gr)
#'
#' @import purrr
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges width
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges Views
#' @importFrom IRanges viewMeans
#' @import GenomeInfoDb
#'
#' @export
coverage_2_matrix <- function(coverage, gr, ncol = 100) {
    stopifnot(is(gr, "GRanges"))
    stopifnot(all(width(gr) >= ncol))
    stopifnot(length(unique(width(gr))) == 1)

    i <- S4Vectors::queryHits(GenomicRanges::findOverlaps(coverage, gr))
    coverage <- coverage[i]
    coverage <- GenomicRanges::coverage(coverage, weight = coverage$score)

    gr <- split(gr, as.character(seqnames(gr)))

    extract_scores <- function(n) {
        binned_gr <- intoNbins(gr[[n]], n = ncol)
        cov <- coverage[[n]]
        view <- Views(cov, start(binned_gr), end(binned_gr))
        matrix(viewMeans(view), ncol = ncol, byrow = TRUE)
    }
    m <- map(names(gr), extract_scores)
    m <- do.call("rbind", m)
    m[is.nan(m)] <- 0
    m
}
