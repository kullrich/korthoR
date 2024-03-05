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
//' @title rcpp_count_mers
//' @name rcpp_count_mers
//' @description returns count all mers
//' @return list
//' @param aavector StringVector [mandatory]
//' @param ncores number of cores [default: 1]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' h <- hiv |>
//'     MSA2dist::subString(s=1, e=12) |>
//'     MSA2dist::cds2aa() |>
//'     as.character()
//' l <- korthoR::rcpp_count_mers(
//'     aavector=h,
//'     ncores=1)
//' l
//' @export rcpp_count_mers
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::List rcpp_count_mers(
  Rcpp::StringVector aavector,
  int ncores=1) {
  int nseq = aavector.size();
  std::vector<std::string> seqnames = aavector.attr("names");
  std::vector<std::string> seq = convertRcppStringVectorToStdVector(aavector);
  std::vector<std::vector<std::string>> out_mers(nseq);
  std::vector<std::vector<int>> out_mers_counts(nseq);
  RcppThread::ProgressBar bar(nseq, 1);
  RcppThread::parallelFor(0, nseq, [&] (int seq_idx) {
    int nsites_seq = aavector[seq_idx].size();
    std::vector<std::string> l_seq;
    std::map<std::string, int> counts_seq;
    for (int j = 0; j < nsites_seq-1; j++) {
      for (int i = 1; i < nsites_seq-j+1; i++) {
        std::string mer;
        mer = seq[seq_idx];
        l_seq.push_back(mer.substr(i-1, j+1));
      }
    }
    for (const auto &s : l_seq) {
      counts_seq[s]++;
    }
    std::vector<std::string> result_seq_mers(counts_seq.size());
    std::vector<int> result_seq_mers_counts(counts_seq.size());
    int result_i = 0;
    for (const auto &pair : counts_seq) {
      result_seq_mers[result_i] = pair.first;
      result_seq_mers_counts[result_i] = pair.second;
      result_i++;
    };
    std::tuple<std::vector<std::string>, std::vector<int>> result_seq;
    result_seq = sort_string_vec_and_int_vec(result_seq_mers, result_seq_mers_counts);
    std::vector<std::string> result_seq_mers_sorted = std::get<0>(result_seq);
    std::vector<int> result_seq_mers_counts_sorted = std::get<1>(result_seq);
    out_mers[seq_idx] = result_seq_mers_sorted;
    out_mers_counts[seq_idx] = result_seq_mers_counts_sorted;
    bar++;
  }, ncores);
  Rcpp::List out(nseq);
  out = iterate_and_create_named_int_vec(out_mers_counts, out_mers, seqnames);
  return out;
}
