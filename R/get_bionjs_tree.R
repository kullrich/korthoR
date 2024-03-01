#' @title get_bionjs_tree
#' @name get_bionjs_tree
#' @description This function take the results from \code{get_jaccard_a_b}
#' and use the distances to calculate a \code{bionjs} tree with the
#' \code{ape} package
#' @param jaccard_df data.frame from \code{get_jaccard_a_b}
#' @param value select value from jaccard result to be used as distances to
#' calculate the bionjs tree [default: jaccard]; choices are "jaccard" or "mash"
#' or "sumdist"
#' @return \code{data.frame}
#' @importFrom stats as.dist
#' @importFrom reshape2 melt dcast
#' @importFrom dplyr filter
#' @importFrom ape bionjs
#' @examples
#' data(hiv, package="MSA2dist")
#' l <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_kmers(k=6)
#' d <- korthoR::get_jaccard_a_b(
#'     kmer_counts_q=l,
#'     kmer_counts_t=l,
#'     k=6)
#' t1 <- korthoR::get_bionjs_tree(jaccard_df=d, value="jaccard")
#' t2 <- korthoR::get_bionjs_tree(jaccard_df=d, value="mash")
#' t3 <- korthoR::get_bionjs_tree(jaccard_df=d, value="sumdist")
#' par(mfrow=c(1,3))
#' plot(t1, main="jaccard")
#' plot(t2, main="mash")
#' plot(t3, main="sumdist")
#' @export get_bionjs_tree
#' @author Kristian K Ullrich

get_bionjs_tree <- function(
    jaccard_df,
    value="jaccard") {
    if(value=="jaccard"){
        j_melted_df <- reshape2::melt(jaccard_df,
        id.vars = c("qname", "tname")) |>
        dplyr::filter(variable=="q1t2_jaccard") |>
        reshape2::dcast(qname~tname, value.var="value")
    } else if(value=="sumdist"){
        j_melted_df <- reshape2::melt(jaccard_df,
        id.vars = c("qname", "tname")) |>
        dplyr::filter(variable=="q1t2_sumdist") |>
        reshape2::dcast(qname~tname, value.var="value")
    } else if(value=="mash"){
        j_melted_df <- reshape2::melt(jaccard_df,
        id.vars = c("qname", "tname")) |>
        dplyr::filter(variable=="q1t2_mash") |>
        reshape2::dcast(qname~tname, value.var="value")
    }
    rownames(j_melted_df) <- j_melted_df[,1]
    j_melted_df <- j_melted_df[,-1]
    if(value=="jaccard"){
        t <- ape::bionjs(as.dist(1-j_melted_df))
    } else if(value=="sumdist"){
        t <- ape::bionjs(as.dist(j_melted_df))
    } else if(value=="mash"){
        t <- ape::bionjs(as.dist(j_melted_df))
    }
    return(t)
}
