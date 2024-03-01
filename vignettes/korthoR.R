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
    korthoR::count_kmers(k=6)
l

## -----------------------------------------------------------------------------
l <- hiv |>
    MSA2dist::cds2aa() |>
    korthoR::count_kmers(k=6, threads=2)
l

