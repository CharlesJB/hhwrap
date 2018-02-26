#' Get narrowPeak filename for demo
#' 
#' @return The complete path to a demo narrowPeak file
#' 
#' @examples
#' peak_file <- get_demo_peaks_file()
#' 
#' @export
get_demo_peaks_file <- function() {
    system.file("extdata/peaks.narrowPeak", package="hmwrap")
}

#' Get bedGraph filenames for demo
#' 
#' @return A vector with the bedGraph filenames used for demo
#' 
#' @examples
#' bg_file <- get_demo_bedGraph_files()
#' 
#' @export
get_demo_bedGraph_files <- function() {
    c(system.file("extdata/file1.bg", package="hmwrap"),
      system.file("extdata/file1.bg", package="hmwrap"))
}
