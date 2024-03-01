#' @title count_mers
#' @name count_mers
#' @description This function count all mers from a character vector or an
#' AAString or an AAStringSet
#' @param aa character vector/AAString/AAStringSet
#' @param threads number of parallel threads [default: 1]
#' @return \code{table} or \code{list} of \code{table}
#' @importFrom utils data
#' @importFrom stringr word
#' @importFrom stats setNames
#' @seealso \code{\link[Biostrings]{XStringSet-class}}
#' @examples
#' data(hiv, package="MSA2dist")
#' l <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_mers()
#' l
#' @export count_mers
#' @author Kristian K Ullrich

count_mers <- function(
    aa,
    threads=1) {
    stopifnot("Error: input needs to be a character or AAString or AAStringSet"=
        any(methods::is(aa, "character") ||
        methods::is(aa, "AAString") || methods::is(aa, "AAStringSet")))
    if(methods::is(aa, "character")){
        if(is.null(names(aa))){
          names(aa) <- "aa"
        }
        korthoR::rcpp_count_mers(aavector=aa, ncores=threads)
    } else if(methods::is(aa, "AAString")){
        korthoR::rcpp_count_mers(aavector=stats::setNames(as.character(aa),
            "aa"), ncores=threads)
    } else {
        if(is.null(names(aa))){
            names(aa) <- seq(from=1, to=length(aa))
        } else {
            names(aa) <- stringr::word(names(aa), 1)
        }
        korthoR::rcpp_count_mers(aavector=as.character(aa), ncores=threads)
    }
}
