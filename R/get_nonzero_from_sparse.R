#' @title get_nonzero_from_sparse
#' @name get_nonzero_from_sparse
#' @description This function count all mers from a character vector
#' @param mm Matrix
#' @return \code{Matrix}
#' @examples
#' mm <- Matrix::Matrix(0, ncol=3, nrow=3, sparse = TRUE)
#' mm[1,1] <- 1
#' mm |> get_nonzero_from_sparse()
#' @export get_nonzero_from_sparse
#' @author Kristian K Ullrich

get_nonzero_from_sparse <- function(
    mm) {
    stopifnot("Error: input needs to be a character"=
        methods::is(mm, "Matrix"))
    non_zero_indices <- which(mm[] != 0, arr.ind = TRUE)
    row_indices <- non_zero_indices[, 1]
    col_indices <- non_zero_indices[, 2]
    values <- mm[non_zero_indices]
    result <- cbind(i=row_indices, j=col_indices, v=values)
    return(result)
}
