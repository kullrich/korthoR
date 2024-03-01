#' @title string2mers
#' @name string2mers
#' @description This function count all mers from a character vector
#' @param string character vector
#' @return \code{table}
#' @examples
#' korthoR::string2mers(
#'     "AACGTGCA")
#' @export string2mers
#' @author Kristian K Ullrich

string2mers <- function(
    string) {
    stopifnot("Error: input needs to be a character"=
        methods::is(string, "character"))
    l <- list()
    n <- nchar(string)
    counter <- 0
    for (j in seq(from=0, to=n-2)) {
        for (i in seq(from=1, to=n-j)) {
            counter <- counter+1
            l[counter] <- substr(string, i, i+j)
        }
    }
    t <- table(unlist(l))
    return(t[order(names(t))])
}
