#' @title split_list_into_chunks
#' @name split_list_into_chunks
#' @description This function splits a list into chunks into at most
#' num_per_chunk elements
#' @param l list to split
#' @param num_per_chunk specify max number per chunk [default: 1000]
#' @return \code{table}
#' @examples
#' korthoR::split_list_into_chunks(
#'     l=vector("list", length=10),
#'     num_per_chunk=2)
#' @export split_list_into_chunks
#' @author Kristian K Ullrich

split_list_into_chunks <- function(
    l,
    num_per_chunk=1000) {
    stopifnot("Error: input needs to be a list"=
        methods::is(l, "list"))
    num_chunks <- ceiling(length(l) / num_per_chunk)
    l_chunks <- vector("list", length=num_chunks)
    for (i in seq_along(l_chunks)) {
        start_index <- (i-1) * num_per_chunk + 1
        end_index <- min(i*num_per_chunk, length(l))
        l_chunks[[i]] <- l[seq(from=start_index, to=end_index)]
    }
    return(l_chunks)
}
