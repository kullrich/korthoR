#' @title count_mers2
#' @name count_mers2
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
#'     korthoR::count_mers2()
#' l
#' @export count_mers2
#' @author Kristian K Ullrich

count_mers2 <- function(
    aa,
    threads=1) {
    stopifnot("Error: input needs to be a character or AAString or AAStringSet"=
        any(methods::is(aa, "character") ||
        methods::is(aa, "AAString") || methods::is(aa, "AAStringSet")))
    if(methods::is(aa, "character")){
        list(string2mers(aa))
    } else if(methods::is(aa, "AAString")){
        list(string2mers(as.character(aa)))
    } else {
        if(is.null(names(aa))){
            names(aa) <- seq(from=1, to=length(aa))
        } else {
            names(aa) <- stringr::word(names(aa), 1)
        }
        lapply(aa, function(x) string2mers(as.character(x)))
    }
}
