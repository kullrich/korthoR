#' @title get_jaccard_a_b
#' @name get_jaccard_a_b
#' @description This function calculates jaccard distance for all pairwise
#' comparison between two kmer_counts list objects.
#' If the parameter k is set, mash distance will be reported.
#' Only jaccard distances will be reported if threshold is met.
#' @param kmer_counts_q kmer_counts
#' @param kmer_counts_t kmer_counts
#' @param k kmer length [default: 6]
#' @param min_jaccard min jaccard distance to report [default: 0.01]
#' @param use_rcpp use rcpp jaccard distance calculation[ default: TRUE]
#' @param use_sparse use rcpp sparse approach [default: TRUE]
#' @param threads number of parallel threads [default: 1]
#' @param num_per_chunk_q specify max number per chunk q [default: 100]
#' @param num_per_chunk_t specify max number per chunk t [default: 1000]
#' @param aa2int convert aa to Int64 [default: FALSE]
#' @return \code{data.frame}
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeForkCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @importFrom stringr word str_c
#' @importFrom dplyr arrange
#' @examples
#' data(hiv, package="MSA2dist")
#' l <- hiv |>
#'     MSA2dist::cds2aa() |>
#'     korthoR::count_kmers(k=6)
#' df <- korthoR::get_jaccard_a_b(
#'     kmer_counts_q=l,
#'     kmer_counts_t=l,
#'     k=6,
#'     min_jaccard=0.01)
#' df
#' # to use multiple threads
#' #df <- korthoR::get_jaccard_a_b(
#' #    kmer_counts_q=l,
#' #    kmer_counts_t=l,
#' #    k=6,
#' #    min_jaccard=0.01,
#' #    threads=2)
#' #df
#' # to use multiple threads and change chunk size
#' #df <- korthoR::get_jaccard_a_b(
#' #    kmer_counts_q=l,
#' #    kmer_counts_t=l,
#' #    k=6,
#' #    min_jaccard=0.01,
#' #    threads=2,
#' #    num_per_chunk=5)
#' #df
#' # to use Int64 kmers
#' # lint64 <- hiv |>
#' #     MSA2dist::cds2aa() |>
#' #     korthoR::count_kmers(
#' #     k=6,
#' #     aa2int=TRUE)
#' #df <- korthoR::get_jaccard_a_b(
#' #    kmer_counts_q=lint64,
#' #    kmer_counts_t=lint64,
#' #    k=6,
#' #    min_jaccard=0.01,
#' #    threads=2,
#' #    num_per_chunk=5,
#' #    aa2int=TRUE)
#' #df
#' @export get_jaccard_a_b
#' @author Kristian K Ullrich

get_jaccard_a_b <- function(
    kmer_counts_q,
    kmer_counts_t,
    k=6,
    min_jaccard=0.01,
    use_rcpp=TRUE,
    use_sparse=TRUE,
    threads=1,
    num_per_chunk_q=100,
    num_per_chunk_t=1000,
    aa2int=FALSE) {
    stopifnot("Error: k needs to be in a range of min 2 to max 12"=
        k %in% seq(from=2, to=12))
    if(!use_rcpp){
        if(.Platform$OS.type == "windows"){
            cl <- parallel::makeCluster(threads)
        }
        if(.Platform$OS.type != "windows"){
            cl <- parallel::makeForkCluster(threads)
        }
        doParallel::registerDoParallel(cl)
        i <- NULL
        j <- NULL
        OUT <- NULL
        OUT <- foreach(i=seq(from=1, to=length(kmer_counts_q)),
            .packages = c('foreach'), .combine=rbind) %dopar% {
            foreach(j=seq(from=1, to=length(kmer_counts_t)),
            .combine=rbind) %do% {
                qname <- setNames(NA, "qname")
                tname <- setNames(NA, "tname")
                q1t2_idx <- NULL
                q2t1_idx <- NULL
                q1t2_uniq <- NULL
                q2t1_uniq <- NULL
                q1t2_idx_sum <- setNames(0, "q1t2_idx_sum")
                q2t1_idx_sum <- setNames(0, "q2t1_idx_sum")
                q1t2_uniq_sum <- setNames(0, "q1t2_uniq_sum")
                q2t1_uniq_sum <- setNames(0, "q2t1_uniq_sum")
                q1t2_jaccard <- setNames(NA, "jaccard")
                q1t2_sumdist <- setNames(NA, "sumdist")
                q1t2_mash <- setNames(NA, "mash")
                q1t2_ani <- setNames(NA, "ani")
                qname <- setNames(names(kmer_counts_q[i]), "qname")
                tname <- setNames(names(kmer_counts_t[j]), "tname")
                q1t2_idx <- which(names(
                    kmer_counts_q[[i]])%in%names(kmer_counts_t[[j]]))
                q2t1_idx <- which(names(
                    kmer_counts_t[[j]])%in%names(kmer_counts_q[[i]]))
                q1t2_uniq <- which(!names(
                    kmer_counts_q[[i]])%in%names(kmer_counts_t[[j]]))
                q2t1_uniq <- which(!names(
                    kmer_counts_t[[j]])%in%names(kmer_counts_q[[i]]))
                q1t2_idx_sum <- setNames(
                    sum(kmer_counts_q[[i]][q1t2_idx]), "q1t2_idx_sum")
                q2t1_idx_sum <- setNames(
                    sum(kmer_counts_t[[j]][q2t1_idx]), "q2t1_idx_sum")
                q1t2_uniq_sum <- setNames(
                    sum(kmer_counts_q[[i]][q1t2_uniq]), "q1t2_uniq_sum")
                q2t1_uniq_sum <- setNames(
                    sum(kmer_counts_t[[j]][q2t1_uniq]), "q2t1_uniq_sum")
                q1t2_jaccard <- setNames(
                    length(q1t2_idx) /
                    (length(q1t2_idx) + length(q1t2_uniq) + length(q2t1_uniq)),
                    "jaccard")
                q1t2_sumdist <- setNames(
                    1-( (q1t2_idx_sum+q2t1_idx_sum) /
                    (q1t2_idx_sum+q2t1_idx_sum+q1t2_uniq_sum+q2t1_uniq_sum) ),
                    "sumdist")
                if(!is.null(k)){
                    q1t2_mash <- setNames(
                    -(1/k) * log( (2*q1t2_jaccard) / (1+q1t2_jaccard) ),
                    "mash")
                    q1t2_ani <- setNames(1-q1t2_mash, "ani")
                }
                if(q1t2_jaccard>min_jaccard){
                    d <- cbind(
                        qname, tname,
                        q1t2_jaccard,
                        q1t2_mash,
                        q1t2_ani,
                        q1t2_sumdist)
                    d
                }
            }
        }
        parallel::stopCluster(cl)
        OUT <- as.data.frame(OUT)
        OUT["q1t2_jaccard"] <- as.numeric(unlist(OUT["q1t2_jaccard"]))
        OUT["q1t2_mash"] <- as.numeric(unlist(OUT["q1t2_mash"]))
        OUT["q1t2_ani"] <- as.numeric(unlist(OUT["q1t2_ani"]))
        OUT["q1t2_sumdist"] <- as.numeric(unlist(OUT["q1t2_sumdist"]))
        rownames(OUT) <- NULL
        return(OUT)
    } else {
        names(kmer_counts_q) <- stringr::str_c("q_",
            stringr::word(names(kmer_counts_q)))
        names(kmer_counts_t) <- stringr::str_c("t_",
            stringr::word(names(kmer_counts_t)))
        kmer_counts_q_chunks <- korthoR::split_list_into_chunks(
            l=kmer_counts_q, num_per_chunk=num_per_chunk_q)
        kmer_counts_t_chunks <- korthoR::split_list_into_chunks(
            l=kmer_counts_t, num_per_chunk=num_per_chunk_t)
        OUT <- NULL
        for(q_chunk in kmer_counts_q_chunks){
            for(t_chunk in kmer_counts_t_chunks){
                if(!use_sparse){
                  if (aa2int) {
                      q_chunk_t_chunk_OUT <- korthoR::rcpp_jaccard_a_b_Int64(
                          kmer_counts_q_Int64=q_chunk,
                          kmer_counts_t_Int64=t_chunk,
                          k=k,
                          min_jaccard=min_jaccard,
                          ncores=threads)
                  } else{
                      q_chunk_t_chunk_OUT <- korthoR::rcpp_jaccard_a_b(
                          kmer_counts_q=q_chunk,
                          kmer_counts_t=t_chunk,
                          k=k,
                          min_jaccard=min_jaccard,
                          ncores=threads)
                  }
                } else {
                    if (aa2int) {
                        q_chunk_t_chunk_OUT <-
                            korthoR::rcpp_jaccard_sparse_a_b_Int64(
                            kmer_counts_q_Int64=q_chunk,
                            kmer_counts_t_Int64=t_chunk,
                            k=k,
                            min_jaccard=min_jaccard,
                            ncores=threads,
                            debug=FALSE)
                    } else {
                        q_chunk_t_chunk_OUT <- korthoR::rcpp_jaccard_sparse_a_b(
                            kmer_counts_q=q_chunk,
                            kmer_counts_t=t_chunk,
                            k=k,
                            min_jaccard=min_jaccard,
                            ncores=threads,
                            debug=FALSE)
                    }
                }
                OUT <- rbind(OUT, q_chunk_t_chunk_OUT)
            }
        }
        OUT[["qname"]] <- sub("q_", "",OUT[["qname"]])
        OUT[["tname"]] <- sub("t_", "", OUT[["tname"]])
        OUT <- dplyr::arrange(OUT, qname, tname)
        OUT <- dplyr::arrange(OUT, qname, tname)
        return(OUT)
    }
}
