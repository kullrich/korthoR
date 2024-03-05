#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
#include <Rcpp.h>
#include <RcppThread.h>
#include <string.h>
#include "kmer_utils.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib korthoR, .registration = TRUE
//' @import Rcpp
//' @title rcpp_count_kmers
//' @name rcpp_count_kmers
//' @description returns count kmers
//' @return list
//' @param aavector StringVector [mandatory]
//' @param k kmer length [default: 6]
//' @param ncores number of cores [default: 1]
//' @param aa2int convert aa to Int64 [default: FALSE]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' h <- hiv |> MSA2dist::cds2aa() |> as.character()
//' l <- korthoR::rcpp_count_kmers(
//'     aavector=h,
//'     k=6,
//'     ncores=1)
//' l
//' # to convert aa to Int64
//' lint64 <- korthoR::rcpp_count_kmers(
//'     aavector=h,
//'     k=6,
//'     ncores=1,
//'     aa2int=TRUE)
//' lint64
//' @export rcpp_count_kmers
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::List rcpp_count_kmers(
  Rcpp::StringVector aavector,
  int k=6,
  int ncores=1,
  bool aa2int=false) {
  if (aa2int) {
    int nseq = aavector.size();
    std::vector<std::string> seqnames = aavector.attr("names");
    std::vector<std::string> seq = convertRcppStringVectorToStdVector(aavector);
    std::vector<std::vector<std::int64_t>> out_kmersInt64(nseq);
    std::vector<std::vector<int>> out_kmers_counts(nseq);
    RcppThread::ProgressBar bar(nseq, 1);
    RcppThread::parallelFor(0, nseq, [&] (int seq_idx) {
      int nsites_seq = seq[seq_idx].size();
      std::vector<std::string> l_seq;
      std::map<std::string, int> counts_seq;
      for (int i = 0; i < nsites_seq-k+1; i++) {
        std::string kmer;
        kmer = seq[seq_idx];
        l_seq.push_back(kmer.substr(i, k));
      }
      for (const auto &s : l_seq) {
        counts_seq[s]++;
      }
      std::vector<std::int64_t> result_seq_kmersInt64(counts_seq.size());
      std::vector<int> result_seq_kmers_counts(counts_seq.size());
      int result_i = 0;
      for (const auto &pair : counts_seq) {
        result_seq_kmersInt64[result_i] = aminoAcidSequenceToInt(pair.first);
        result_seq_kmers_counts[result_i] = pair.second;
        result_i++;
      };
      std::tuple<std::vector<std::int64_t>, std::vector<int>> result_seq;
      result_seq = sort_string_vec_and_int_vecInt64(result_seq_kmersInt64, result_seq_kmers_counts);
      std::vector<std::int64_t> result_seq_kmers_sorted = std::get<0>(result_seq);
      std::vector<int> result_seq_kmers_counts_sorted = std::get<1>(result_seq);
      out_kmersInt64[seq_idx] = result_seq_kmers_sorted;
      out_kmers_counts[seq_idx] = result_seq_kmers_counts_sorted;
      bar++;
    }, ncores);
    Rcpp::List out(nseq);
    out = iterate_and_create_named_int_vecInt64(out_kmers_counts, out_kmersInt64, seqnames);
    return out;
  } else {
    int nseq = aavector.size();
    std::vector<std::string> seqnames = aavector.attr("names");
    std::vector<std::string> seq = convertRcppStringVectorToStdVector(aavector);
    std::vector<std::vector<std::string>> out_kmers(nseq);
    std::vector<std::vector<int>> out_kmers_counts(nseq);
    RcppThread::ProgressBar bar(nseq, 1);
    RcppThread::parallelFor(0, nseq, [&] (int seq_idx) {
      int nsites_seq = seq[seq_idx].size();
      std::vector<std::string> l_seq;
      std::map<std::string, int> counts_seq;
      for (int i = 0; i < nsites_seq-k+1; i++) {
        std::string kmer;
        kmer = seq[seq_idx];
        l_seq.push_back(kmer.substr(i, k));
      }
      for (const auto &s : l_seq) {
        counts_seq[s]++;
      }
      std::vector<std::string> result_seq_kmers(counts_seq.size());
      std::vector<int> result_seq_kmers_counts(counts_seq.size());
      int result_i = 0;
      for (const auto &pair : counts_seq) {
        result_seq_kmers[result_i] = pair.first;
        result_seq_kmers_counts[result_i] = pair.second;
        result_i++;
      };
      std::tuple<std::vector<std::string>, std::vector<int>> result_seq;
      result_seq = sort_string_vec_and_int_vec(result_seq_kmers, result_seq_kmers_counts);
      std::vector<std::string> result_seq_kmers_sorted = std::get<0>(result_seq);
      std::vector<int> result_seq_kmers_counts_sorted = std::get<1>(result_seq);
      out_kmers[seq_idx] = result_seq_kmers_sorted;
      out_kmers_counts[seq_idx] = result_seq_kmers_counts_sorted;
      bar++;
    }, ncores);
    Rcpp::List out(nseq);
    out = iterate_and_create_named_int_vec(out_kmers_counts, out_kmers, seqnames);
    return out;
  }
}
