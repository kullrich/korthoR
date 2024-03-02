## ----echo=FALSE, results='hide', warning=FALSE, message=FALSE-----------------
suppressPackageStartupMessages({
    library(korthoR)
    library(MSA2dist)
    library(Biostrings)
    library(ape)
    library(dplyr)
    library(reshape2)
    })

## -----------------------------------------------------------------------------
# load korthoR
library(korthoR)
# load example data
data(hiv, package="MSA2dist")

## -----------------------------------------------------------------------------
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6)
l

## -----------------------------------------------------------------------------
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6,
        threads=2)
l

## -----------------------------------------------------------------------------
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6,
        threads=2,
        aa2int=TRUE)
l

## -----------------------------------------------------------------------------
system.time(d <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(k=6) |>
    korthoR::get_jaccard_from_self(
        k=6))
d

## -----------------------------------------------------------------------------
system.time(d <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6,
        threads=2) |>
    korthoR::get_jaccard_from_self(
        k=6,
        threads=2))
d

## -----------------------------------------------------------------------------
system.time(d <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6,
        threads=2,
        aa2int=TRUE) |>
    korthoR::get_jaccard_from_self(
        k=6,
        threads=2,
        aa2int=TRUE))
d

## -----------------------------------------------------------------------------
system.time(d <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
        k=6,
        threads=2) |>
    korthoR::get_jaccard_from_self(
        k=6,
        use_rcpp = TRUE,
        use_sparse = FALSE,
        threads=2))
d

## -----------------------------------------------------------------------------
system.time(d <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
      k=6,
      threads=2) |>
    korthoR::get_jaccard_from_self(
        k=6,
        use_rcpp = FALSE,
        use_sparse = FALSE,
        threads=2))
d

## -----------------------------------------------------------------------------
data(hiv, package="MSA2dist")
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(
      k=6)
d <- korthoR::get_jaccard_a_b(
    kmer_counts_q=l,
    kmer_counts_t=l,
    k=6)
t1 <- korthoR::get_bionjs_tree(
    jaccard_df=d,
    value="jaccard")
t2 <- korthoR::get_bionjs_tree(
    jaccard_df=d,
    value="mash")
t3 <- korthoR::get_bionjs_tree(
    jaccard_df=d,
    value="sumdist")
par(mfrow=c(1,3))
plot(t1, main="jaccard")
plot(t2, main="mash")
plot(t3, main="sumdist")

## -----------------------------------------------------------------------------
# Download peptides
PAO1 <- Biostrings::readAAStringSet("https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.faa.gz")
SBW25 <- Biostrings::readAAStringSet("https://www.pseudomonas.com/downloads/pseudomonas/pgd_r_22_1/Pseudomonas_fluorescens_SBW25_116/Pseudomonas_fluorescens_SBW25_116.faa.gz")
# Count kmers using multiple threads converting AA to Int64
l_PAO1_int64 <- korthoR::count_kmers(
  aa=PAO1,
  k=6,
  threads=8,
  aa2int=TRUE)
l_SBW25_int64 <- korthoR::count_kmers(
  aa=SBW25,
  k=6,
  threads=8,
  aa2int=TRUE)
# Get jaccard distance using multiple threads and rcpp converting AA to Int64
system.time(df_PAO1_SBW25 <- korthoR::get_jaccard_a_b(
    kmer_counts_q=l_PAO1_int64,
    kmer_counts_t=l_SBW25_int64,
    k=6,
    min_jaccard=0.01,
    threads=8,
    aa2int=TRUE))
head(df_PAO1_SBW25)

## ----sessionInfo, echo=TRUE---------------------------------------------------
sessionInfo()

