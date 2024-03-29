% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_jaccard_a_b}
\alias{rcpp_jaccard_a_b}
\title{rcpp_jaccard_a_b}
\arguments{
\item{kmer_counts_q}{list [mandatory]}

\item{kmer_counts_t}{list [mandatory]}

\item{k}{kmer length [default: 6]}

\item{min_jaccard}{min jaccard distance to report [default: 0.01]}

\item{ncores}{number of cores [default: 1]}
}
\value{
list
}
\description{
returns jaccard
}
\examples{
## load example sequence data
data("hiv", package="MSA2dist")
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(k=6)
d <- korthoR::rcpp_jaccard_a_b(
    kmer_counts_q=l,
    kmer_counts_t=l,
    k=6,
    min_jaccard=0.01,
    ncores=1)
d
}
\author{
Kristian K Ullrich
}
