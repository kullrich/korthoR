#' @title string2kmers
#' @name string2kmers
#' @description This function count kmers from a character
#' @param string character vector
#' @param k kmer length [default: 6]
#' @return \code{table}
#' @examples
#' korthoR::string2kmers(
#'     "AACGTGCA",
#'     k=6)
#' @export string2kmers
#' @author Kristian K Ullrich

string2kmers <- function(
    string,
    k=6) {
    stopifnot("Error: input needs to be a character"=
        methods::is(string, "character"))
    stopifnot("Error: k needs to be in a range of min 2 to max 12"=
        k %in% seq(from=2, to=12))
    l <- list()
    n <- nchar(string)
    for (i in seq(from=1, to=n-k+1)) {
        l[i] <- substr(string, i, i+k-1)
    }
    t <- table(unlist(l))
    return(t[order(names(t))])
}
