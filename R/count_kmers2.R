#' @title count_kmers2
#' @name count_kmers2
#' @description This function count kmers from a character vector or an
#' AAString or an AAStringSet
#' @param aa character vector/AAString/AAStringSet
#' @param k kmer length [default: 6]
#' @return \code{table} or \code{list} of \code{table}
#' @importFrom utils data
#' @importFrom stringr word
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' data(hiv, package="MSA2dist")
#' l <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_kmers2(k=6)
#' l
#' @export count_kmers2
#' @author Kristian K Ullrich

count_kmers2 <- function(
    aa,
    k=6) {
    stopifnot("Error: input needs to be a character or AAString or AAStringSet"=
        any(methods::is(aa, "character") ||
        methods::is(aa, "AAString") || methods::is(aa, "AAStringSet")))
    stopifnot("Error: k needs to be in a range of min 2 to max 12"=
        k %in% seq(from=2, to=12))
    if(methods::is(aa, "character")){
        if(length(aa)>1){
            lapply(aa, function(x) string2kmers(x, k=k))
        } else {
            list(string2kmers(aa, k=k))
        }
    } else if(methods::is(aa, "AAString")){
        list(string2kmers(as.character(aa), k=k))
    } else {
        if(is.null(names(aa))){
            names(aa) <- seq(from=1, to=length(aa))
        } else {
            names(aa) <- stringr::word(names(aa), 1)
        }
        lapply(aa, function(x) string2kmers(as.character(x), k=k))
    }
}
