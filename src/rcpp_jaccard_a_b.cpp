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
//' @title rcpp_jaccard_a_b
//' @name rcpp_jaccard_a_b
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q list [mandatory]
//' @param kmer_counts_t list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param ncores number of cores [default: 1]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' l <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(k=6)
//' d <- korthoR::rcpp_jaccard_a_b(
//'     kmer_counts_q=l,
//'     kmer_counts_t=l,
//'     k=6,
//'     min_jaccard=0.01,
//'     ncores=1)
//' d
//' @export rcpp_jaccard_a_b
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_a_b(
  Rcpp::List kmer_counts_q,
  Rcpp::List kmer_counts_t,
  int k=6,
  double min_jaccard=0.01,
  int ncores=1) {
  int nseq_q = kmer_counts_q.size();
  int nseq_t = kmer_counts_t.size();
  std::vector<std::string> seqnames_q = kmer_counts_q.attr("names");
  std::vector<std::string> seqnames_t = kmer_counts_t.attr("names");
  std::vector<std::vector<std::string>> seq_q_kmers_sorted(nseq_q);
  std::vector<std::vector<std::string>> seq_t_kmers_sorted(nseq_t);
  std::vector<std::vector<int>> seq_q_kmers_counts_sorted(nseq_q);
  std::vector<std::vector<int>> seq_t_kmers_counts_sorted(nseq_t);
  seq_q_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  seq_t_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  seq_q_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  seq_t_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  std::vector<std::vector<std::vector<double>>> seqj_seqk_distances(
    nseq_q,
    std::vector<std::vector<double>>(
      nseq_t,
      std::vector<double>(
        6,
        0.0)));
  RcppThread::ThreadPool pool(ncores);
  pool.parallelFor(0, nseq_q, [&] (int seq_q_i) {
    pool.parallelFor(0, nseq_t, [&, seq_q_i] (int seq_t_i) {
      std::vector<std::string> seq_q_i_kmers_sorted = seq_q_kmers_sorted[seq_q_i];
      std::vector<int> seq_q_i_kmers_counts_sorted = seq_q_kmers_counts_sorted[seq_q_i];
      std::vector<std::string> seq_t_i_kmers_sorted = seq_t_kmers_sorted[seq_t_i];
      std::vector<int> seq_t_i_kmers_counts_sorted = seq_t_kmers_counts_sorted[seq_t_i];
      seqj_seqk_distances[seq_q_i][seq_t_i] = getJaccardByIntegerVector(
        seq_q_i_kmers_sorted,
        seq_t_i_kmers_sorted,
        seq_q_i_kmers_counts_sorted,
        seq_t_i_kmers_counts_sorted,
        seq_q_i, seq_t_i,
        k);
    });
  });
  pool.join();
  std::vector<std::string> out_qname;
  std::vector<std::string> out_tname;
  std::vector<double> out_jaccard;
  std::vector<double> out_mash;
  std::vector<double> out_ani;
  std::vector<double> out_sumdist;
  for (const auto& row : seqj_seqk_distances) {
    for (const auto& distances : row) {
      if ( distances[2] > min_jaccard ) {
        out_qname.push_back(seqnames_q[distances[0]]);
        out_tname.push_back(seqnames_t[distances[1]]);
        out_jaccard.push_back(distances[2]);
        out_mash.push_back(distances[3]);
        out_ani.push_back(distances[4]);
        out_sumdist.push_back(distances[5]);
      }
    }
  }
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
  Rcpp::Named("qname") = out_qname,
  Rcpp::Named("tname") = out_tname,
  Rcpp::Named("q1t2_jaccard") = out_jaccard,
  Rcpp::Named("q1t2_mash") = out_mash,
  Rcpp::Named("q1t2_ani") = out_ani,
  Rcpp::Named("q1t2_sumdist") = out_sumdist);
  return df;
}
