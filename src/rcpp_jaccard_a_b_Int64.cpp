#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <string.h>
#include "kmer_utils.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib korthoR, .registration = TRUE
//' @import Rcpp
//' @title rcpp_jaccard_a_b_Int64
//' @name rcpp_jaccard_a_b_Int64
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q_Int64 list [mandatory]
//' @param kmer_counts_t_Int64 list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param ncores number of cores [default: 1]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' lint64 <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(
//'         k=6,
//'         aa2int=TRUE)
//' d <- korthoR::rcpp_jaccard_a_b_Int64(
//'     kmer_counts_q_Int64=lint64,
//'     kmer_counts_t_Int64=lint64,
//'     k=6,
//'     min_jaccard=0.01,
//'     ncores=1)
//' d
//' @export rcpp_jaccard_a_b_Int64
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_a_b_Int64(
  Rcpp::List kmer_counts_q_Int64,
  Rcpp::List kmer_counts_t_Int64,
  int k=6,
  double min_jaccard=0.01,
  int ncores=1) {
  int nseq_q = kmer_counts_q_Int64.size();
  int nseq_t = kmer_counts_t_Int64.size();
  std::vector<std::string> seqnames_q = kmer_counts_q_Int64.attr("names");
  std::vector<std::string> seqnames_t = kmer_counts_t_Int64.attr("names");
  std::vector<std::vector<std::int64_t>> seq_q_kmers_sorted(nseq_q);
  std::vector<std::vector<std::int64_t>> seq_t_kmers_sorted(nseq_t);
  std::vector<std::vector<int>> seq_q_kmers_counts_sorted(nseq_q);
  std::vector<std::vector<int>> seq_t_kmers_counts_sorted(nseq_t);
  seq_q_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVectorInt64(kmer_counts_q_Int64);
  seq_t_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVectorInt64(kmer_counts_t_Int64);
  seq_q_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_q_Int64);
  seq_t_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_t_Int64);
  std::vector<std::vector<std::vector<double>>> seqj_seqk_distances(
    nseq_q,
    std::vector<std::vector<double>>(
      nseq_t,
      std::vector<double>(
        6,
        0.0)));
  RcppThread::ThreadPool pool(ncores);
  //RcppThread::parallelFor(0, nseq_q, [&] (int seq_q_i) {
  pool.parallelFor(0, nseq_q, [&] (int seq_q_i) {
    //for (int seq_t_i = 0; seq_t_i < nseq_t; ++seq_t_i) {
    pool.parallelFor(0, nseq_t, [&, seq_q_i] (int seq_t_i) {
      std::vector<std::int64_t> seq_q_i_kmers_sorted = seq_q_kmers_sorted[seq_q_i];
      std::vector<int> seq_q_i_kmers_counts_sorted = seq_q_kmers_counts_sorted[seq_q_i];
      std::vector<std::int64_t> seq_t_i_kmers_sorted = seq_t_kmers_sorted[seq_t_i];
      std::vector<int> seq_t_i_kmers_counts_sorted = seq_t_kmers_counts_sorted[seq_t_i];
      seqj_seqk_distances[seq_q_i][seq_t_i] = getJaccardByIntegerVectorInt64(
        seq_q_i_kmers_sorted,
        seq_t_i_kmers_sorted,
        seq_q_i_kmers_counts_sorted,
        seq_t_i_kmers_counts_sorted,
        seq_q_i, seq_t_i,
        k);
    });
  });
  pool.join();
      //}
  //}, ncores);
  std::vector<std::string> seqjnamesVector;
  std::vector<std::string> seqknamesVector;
  std::vector<double> seqj_seqk_jaccardVector;
  std::vector<double> seqj_seqk_mashVector;
  std::vector<double> seqj_seqk_aniVector;
  std::vector<double> seqj_seqk_sumdistVector;
  for (const auto& row : seqj_seqk_distances) {
    for (const auto& distances : row) {
      if ( distances[2] > min_jaccard ) {
        seqjnamesVector.push_back(seqnames_q[distances[0]]);
        seqknamesVector.push_back(seqnames_t[distances[1]]);
        seqj_seqk_jaccardVector.push_back(distances[2]);
        seqj_seqk_mashVector.push_back(distances[3]);
        seqj_seqk_aniVector.push_back(distances[4]);
        seqj_seqk_sumdistVector.push_back(distances[5]);
      }
    }
  }
  Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("qname") = seqjnamesVector,
    Named("tname") = seqknamesVector,
    Named("q1t2_jaccard") = seqj_seqk_jaccardVector,
    Named("q1t2_mash") = seqj_seqk_mashVector,
    Named("q1t2_ani") = seqj_seqk_aniVector,
    Named("q1t2_sumdist") = seqj_seqk_sumdistVector);
  return df;
}
