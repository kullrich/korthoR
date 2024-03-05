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
//' @title rcpp_jaccard_sparse_a_b
//' @name rcpp_jaccard_sparse_a_b
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q list [mandatory]
//' @param kmer_counts_t list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param ncores number of cores [default: 1]
//' @param debug print to console [default: FALSE]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' l <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(k=6)
//' d <- korthoR::rcpp_jaccard_sparse_a_b(
//'     kmer_counts_q=l,
//'     kmer_counts_t=l,
//'     k=6,
//'     min_jaccard=0.01,
//'     ncores=1,
//'     debug=FALSE)
//' d
//' @export rcpp_jaccard_sparse_a_b
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_sparse_a_b(
  Rcpp::List kmer_counts_q,
  Rcpp::List kmer_counts_t,
  int k=6,
  double min_jaccard=0.01,
  int ncores=1,
  bool debug=false) {
  int nseq_q = kmer_counts_q.size();
  int nseq_t = kmer_counts_t.size();
  //create std::vector for names, kmer counts, kmers
  //names
  std::vector<std::string> seqnames_q = kmer_counts_q.attr("names");
  std::vector<std::string> seqnames_t = kmer_counts_t.attr("names");
  //kmers
  std::vector<std::vector<std::string>> seq_q_kmers_sorted(nseq_q);
  std::vector<std::vector<std::string>> seq_t_kmers_sorted(nseq_t);
  //kmer counts
  std::vector<std::vector<int>> seq_q_kmers_counts_sorted(nseq_q);
  std::vector<std::vector<int>> seq_t_kmers_counts_sorted(nseq_t);
  //convert Rcpp vectors into std::vector
  seq_q_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  seq_t_kmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  seq_q_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  seq_t_kmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  //create kmerMaps
  std::map<std::string, std::vector<int>> kmerMap_q;
  std::map<std::string, std::vector<int>> kmerMap_t;
  //create a 2D boolean array to check if a comparison has already been done
  std::vector<std::vector<bool>> comparisonResults(nseq_q, std::vector<bool>(nseq_t, false));
  //create vector to keep separate maps
  auto start_QkmerMap = std::chrono::steady_clock::now();
  std::vector<std::map<std::string, int>> kmerMap_qs(nseq_q);
  kmerMap_qs.resize(nseq_q);
  RcppThread::parallelFor(0, nseq_q, [&] (int q_i) {
    //create individual kmerMap_q_i
    std::map<std::string, int> kmerMap_q_i;
    //create temporary seq_q_i_kmers_sorted
    std::vector<std::string> seq_q_i_kmers_sorted = seq_q_kmers_sorted[q_i];
    //convert seq_q_i_kmers_sorted into map
    kmerMap_q_i = vectorToMap(seq_q_i_kmers_sorted, q_i);
    kmerMap_qs[q_i] = kmerMap_q_i;
  }, ncores);
  //combine maps
  kmerMap_q = combineMaps(kmerMap_qs);
  auto end_QkmerMap = std::chrono::steady_clock::now();
  auto duration_QkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_QkmerMap - start_QkmerMap);
  if (debug) {
    std::cout << "Time taken: QkmerMap creation " << duration_QkmerMap.count() << " milliseconds" << std::endl;
    std::cout << "QkmerMap_n size: " << kmerMap_qs.size() << std::endl;
    std::cout << "number of kmers QkmerMap: " << kmerMap_q.size() << std::endl;
  }
  //create vector to keep separate maps
  auto start_TkmerMap = std::chrono::steady_clock::now();
  std::vector<std::map<std::string, int>> kmerMap_ts(nseq_t);
  kmerMap_ts.resize(nseq_t);
  RcppThread::parallelFor(0, nseq_t, [&] (int t_i) {
    //create individual kmerMap_t_i
    std::map<std::string, int> kmerMap_t_i;
    //create temporary seq_t_i_kmers_sorted
    std::vector<std::string> seq_t_i_kmers_sorted = seq_t_kmers_sorted[t_i];
    //convert seq_t_i_kmers_sorted into map
    kmerMap_t_i = vectorToMap(seq_t_i_kmers_sorted, t_i);
    kmerMap_ts[t_i] = kmerMap_t_i;
  }, ncores);
  //combine maps
  kmerMap_t = combineMaps(kmerMap_ts);
  auto end_TkmerMap = std::chrono::steady_clock::now();
  auto duration_TkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_TkmerMap - start_TkmerMap);
  if (debug) {
    std::cout << "Time taken: TkmerMap creation " << duration_TkmerMap.count() << " milliseconds" << std::endl;
    std::cout << "TkmerMap_n size: " << kmerMap_ts.size() << std::endl;
    std::cout << "number of kmers TkmerMap: " << kmerMap_t.size() << std::endl;
  }
  //direct jaccard
  auto start_calcDist = std::chrono::steady_clock::now();
  std::vector<std::vector<double>> jaccard_distances;
  jaccard_distances = directJaccard(
    kmerMap_q,
    kmerMap_t,
    comparisonResults,
    seq_q_kmers_sorted,
    seq_t_kmers_sorted,
    seq_q_kmers_counts_sorted,
    seq_t_kmers_counts_sorted,
    k,
    min_jaccard);
  auto end_calcDist = std::chrono::steady_clock::now();
  auto duration_calcDist = std::chrono::duration_cast<std::chrono::milliseconds>(end_calcDist - start_calcDist);
  if (debug) {
    std::cout << "Time taken: distance calculation " << duration_calcDist.count() << " milliseconds" << std::endl;
  }
  //create output vectors
  std::vector<std::string> out_qname;
  std::vector<std::string> out_tname;
  std::vector<double> out_jaccard;
  std::vector<double> out_mash;
  std::vector<double> out_ani;
  std::vector<double> out_sumdist;
  for (const auto& jd_i : jaccard_distances) {
    out_qname.push_back(seqnames_q[static_cast<int>(jd_i[0])]);
    out_tname.push_back(seqnames_t[static_cast<int>(jd_i[1])]);
    out_jaccard.push_back(jd_i[2]);
    out_mash.push_back(jd_i[3]);
    out_ani.push_back(jd_i[4]);
    out_sumdist.push_back(jd_i[5]);
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
