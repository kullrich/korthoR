#' @title count_kmers
#' @name count_kmers
#' @description This function count kmers from a character vector or an
#' AAString or an AAStringSet
#' @param aa character vector/AAString/AAStringSet
#' @param k kmer length [default: 6]
#' @param threads number of parallel threads [default: 1]
#' @param aa2int convert aa to Int64 [default: FALSE]
#' @return \code{table} or \code{list} of \code{table}
#' @importFrom utils data
#' @importFrom stringr word
#' @importFrom stats setNames
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' data(hiv, package="MSA2dist")
#' l <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_kmers(k=6)
#' l
#' # to convert aa to Int64
#' lint64 <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_kmers(
#'         k=6,
#'         aa2int=TRUE)
#' lint64
#' @export count_kmers
#' @author Kristian K Ullrich

count_kmers <- function(
    aa,
    k=6,
    threads=1,
    aa2int=FALSE) {
    stopifnot("Error: input needs to be a character or AAString or AAStringSet"=
        any(methods::is(aa, "character") ||
        methods::is(aa, "AAString") || methods::is(aa, "AAStringSet")))
    stopifnot("Error: k needs to be in a range of min 2 to max 12"=
        k %in% seq(from=2, to=12))
    if(methods::is(aa, "character")){
        if(is.null(names(aa))){
          names(aa) <- paste0("aa", seq(from=1, to=length(aa)))
        }
        korthoR::rcpp_count_kmers(
            aavector=aa,
            k=k,
            ncores=threads,
            aa2int=aa2int)
    } else if(methods::is(aa, "AAString")){
        korthoR::rcpp_count_kmers(
            aavector=stats::setNames(as.character(aa), "aa"),
            k=k,
            ncores=threads,
            aa2int=aa2int)
    } else {
        if(is.null(names(aa))){
            names(aa) <- seq(from=1, to=length(aa))
        } else {
            names(aa) <- stringr::word(names(aa), 1)
        }
        korthoR::rcpp_count_kmers(
            aavector=as.character(aa),
            k=k,
            ncores=threads,
            aa2int=aa2int)
    }
}
