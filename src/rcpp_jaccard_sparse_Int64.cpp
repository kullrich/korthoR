//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <string.h>
#include "kmer_utils.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' @useDynLib korthoR, .registration = TRUE
//' @import Rcpp
//' @title rcpp_jaccard_sparse_Int64
//' @name rcpp_jaccard_sparse_Int64
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q_Int64 list [mandatory]
//' @param kmer_counts_t_Int64 list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param ncores number of cores [default: 1]
//' @param debug print to console [default: FALSE]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' l <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(
//'         k=6,
//'         aa2int=TRUE)
//' d <- korthoR::rcpp_jaccard_sparse_Int64(
//'     kmer_counts_q=l,
//'     kmer_counts_t=l,
//'     k=6,
//'     min_jaccard=0.01,
//'     ncores=1,
//'     debug=FALSE)
//' d
//' @export rcpp_jaccard_sparse_Int64
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_sparse_Int64(
  Rcpp::List kmer_counts_q_Int64,
  Rcpp::List kmer_counts_t_Int64,
  int k=6,
  double min_jaccard=0.01,
  int ncores=1,
  bool debug=false) {
  //get sizes
  int nseq_q = kmer_counts_q_Int64.size();
  int nseq_t = kmer_counts_t_Int64.size();
  //create std::vector for names, kmer counts, kmers
  //names
  std::vector<std::string> seqnames_q = kmer_counts_q_Int64.attr("names");
  std::vector<std::string> seqnames_t = kmer_counts_t_Int64.attr("names");
  //kmers
  std::vector<std::vector<std::int64_t>> seq_q_kmers_sorted(nseq_q);
  std::vector<std::vector<std::int64_t>> seq_t_kmers_sorted(nseq_t);
  //kmer counts
  std::vector<std::vector<int>> seq_q_kmers_counts_sorted(nseq_q);
  std::vector<std::vector<int>> seq_t_kmers_counts_sorted(nseq_t);
  //convert Rcpp vectors into std::vector
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
  //create kmerMaps
  std::map<std::int64_t, std::vector<int>> kmerMap_q;
  std::map<std::int64_t, std::vector<int>> kmerMap_t;
  std::map<std::pair<int, int>, bool> candidatePairs;
  //create vector to keep separate maps
  std::vector<std::map<std::int64_t, int>> kmerMap_qs(nseq_q);
  kmerMap_qs.resize(nseq_q);
  RcppThread::parallelFor(0, nseq_q, [&] (int q_i) {
    std::map<std::int64_t, int> kmerMap_q_i;
    std::vector<std::int64_t> seq_q_i_kmers_sorted = seq_q_kmers_sorted[q_i];
    kmerMap_q_i = vectorToInt64Map(seq_q_i_kmers_sorted, q_i);
    kmerMap_qs[q_i] = kmerMap_q_i;
  }, ncores);
  //combine maps
  kmerMap_q = combineInt64Maps(kmerMap_qs);
  //create vector to keep separate maps
  std::vector<std::map<std::int64_t, int>> kmerMap_ts(nseq_t);
  kmerMap_ts.resize(nseq_t);
  RcppThread::parallelFor(0, nseq_t, [&] (int t_i) {
    std::map<std::int64_t, int> kmerMap_t_i;
    std::vector<std::int64_t> seq_t_i_kmers_sorted = seq_t_kmers_sorted[t_i];
    kmerMap_t_i = vectorToInt64Map(seq_t_i_kmers_sorted, t_i);
    kmerMap_ts[t_i] = kmerMap_t_i;
  }, ncores);
  //combine maps
  kmerMap_t = combineInt64Maps(kmerMap_ts);
  if (kmerMap_q.size() < kmerMap_t.size()) {
    for (const auto& kmerQ : kmerMap_q) {
    //RcppThread::parallelFor(0, kmerMap_q.size(), [&](int i) {
      //auto it = std::next(std::begin(kmerMap_q), i);
      //const auto& kmerQ = *it;
      auto kmerT_it = kmerMap_t.lower_bound(kmerQ.first);
      if (kmerT_it != kmerMap_t.end() && !(kmerMap_t.key_comp()(kmerQ.first, kmerT_it->first))) {
        for (const auto& seq_q_i : kmerQ.second) {
          for (const auto& seq_t_i : kmerT_it->second) {
            std::pair<int, int> seq_q_seq_t_pair(seq_q_i, seq_t_i);
            auto it_pair = candidatePairs.find(seq_q_seq_t_pair);
            if (it_pair == candidatePairs.end()) {
              seqj_seqk_distances[seq_q_i][seq_t_i] = getJaccardByIntegerVectorInt64(
                seq_q_kmers_sorted[seq_q_i],
                seq_t_kmers_sorted[seq_t_i],
                seq_q_kmers_counts_sorted[seq_q_i],
                seq_t_kmers_counts_sorted[seq_t_i],
                seq_q_i,
                seq_t_i,
                k);
              candidatePairs[seq_q_seq_t_pair] = true;
            }
          }
        }
      }
    //}, ncores);
    }
  }  else {
    for (const auto& kmerT : kmerMap_t) {
    //RcppThread::parallelFor(0, kmerMap_t.size(), [&](int j) {
      //auto it = std::next(std::begin(kmerMap_t), j);
      //const auto& kmerT = *it;
      auto kmerQ_it = kmerMap_q.lower_bound(kmerT.first);
      if (kmerQ_it != kmerMap_q.end() && !(kmerMap_q.key_comp()(kmerT.first, kmerQ_it->first))) {
        for (const auto& seq_q_i : kmerQ_it->second) {
          for (const auto& seq_t_i : kmerT.second) {
            std::pair<int, int> seq_q_seq_t_pair(seq_q_i, seq_t_i);
            auto it_pair = candidatePairs.find(seq_q_seq_t_pair);
            if (it_pair == candidatePairs.end()) {
              seqj_seqk_distances[seq_q_i][seq_t_i] = getJaccardByIntegerVectorInt64(
                seq_q_kmers_sorted[seq_q_i],
                seq_t_kmers_sorted[seq_t_i],
                seq_q_kmers_counts_sorted[seq_q_i],
                seq_t_kmers_counts_sorted[seq_t_i],
                seq_q_i,
                seq_t_i,
                k);
              candidatePairs[seq_q_seq_t_pair] = true;
            }
          }
        }
      }
    //}, ncores);
    }
  }
  //create output vectors
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
  //create Rcpp::DataFrame
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
  Rcpp::Named("qname") = out_qname,
  Rcpp::Named("tname") = out_tname,
  Rcpp::Named("q1t2_jaccard") = out_jaccard,
  Rcpp::Named("q1t2_mash") = out_mash,
  Rcpp::Named("q1t2_ani") = out_ani,
  Rcpp::Named("q1t2_sumdist") = out_sumdist);
  //Rcpp::DataFrame df;
  return df;
}
