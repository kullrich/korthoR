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
//' @title rcpp_jaccard_sparse_a_b_Int64
//' @name rcpp_jaccard_sparse_a_b_Int64
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q_Int64 list [mandatory]
//' @param kmer_counts_t_Int64_Int64 list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param ncores number of cores [default: 1]
//' @param debug print to console [default: FALSE]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' lint64 <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(
//'         k=6,
//'         aa2int=TRUE)
//' d <- korthoR::rcpp_jaccard_sparse_a_b_Int64(
//'     kmer_counts_q_Int64=lint64,
//'     kmer_counts_t_Int64=lint64,
//'     k=6,
//'     min_jaccard=0.01,
//'     ncores=1,
//'     debug=FALSE)
//' d
//' @export rcpp_jaccard_sparse_a_b_Int64
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_sparse_a_b_Int64(
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
  //create kmerMaps
  std::map<std::int64_t, std::vector<int>> kmerMap_q;
  std::map<std::int64_t, std::vector<int>> kmerMap_t;
  //create a 2D boolean array to check if a comparison has already been done
  std::vector<std::vector<bool>> comparisonResults(nseq_q, std::vector<bool>(nseq_t, false));
  // Define a mutex to protect access to kmerMap_q and kmerMap_t
  //std::mutex mtx_q;
  //std::mutex mtx_t;
  //create vector to keep separate maps
  std::vector<std::map<std::int64_t, int>> kmerMap_qs(nseq_q);
  kmerMap_qs.resize(nseq_q);
  RcppThread::parallelFor(0, nseq_q, [&] (int q_i) {
    //create individual kmerMap_q_i
    std::map<std::int64_t, int> kmerMap_q_i;
    //create temporary seq_q_i_kmers_sorted
    std::vector<std::int64_t> seq_q_i_kmers_sorted = seq_q_kmers_sorted[q_i];
    //convert seq_q_i_kmers_sorted into map
    kmerMap_q_i = vectorToInt64Map(seq_q_i_kmers_sorted, q_i);
    //for (int kmer_i = 0; kmer_i < seq_q_i_kmers_sorted.size(); ++kmer_i) {
      // Protect access to kmerMap_qs with a mutex
      //std::lock_guard<std::mutex> lock(mtx_q);
      //fill kmerMap_q with kmers and link them to sequence q id
      //kmerMap_q[seq_q_i_kmers_sorted[kmer_i]].push_back(q_i);
      //kmerMap_q_i[seq_q_i_kmers_sorted[kmer_i]]=q_i;
    //}
    kmerMap_qs[q_i] = kmerMap_q_i;
  }, ncores);
  //combine maps
  kmerMap_q = combineInt64Maps(kmerMap_qs);
  //create vector to keep separate maps
  std::vector<std::map<std::int64_t, int>> kmerMap_ts(nseq_t);
  kmerMap_ts.resize(nseq_t);
  RcppThread::parallelFor(0, nseq_t, [&] (int t_i) {
    //create individual kmerMap_t_i
    std::map<std::int64_t, int> kmerMap_t_i;
    //create temporary seq_t_i_kmers_sorted
    std::vector<std::int64_t> seq_t_i_kmers_sorted = seq_t_kmers_sorted[t_i];
    //convert seq_t_i_kmers_sorted into map
    kmerMap_t_i = vectorToInt64Map(seq_t_i_kmers_sorted, t_i);
    //for (int kmer_i = 0; kmer_i < seq_t_i_kmers_sorted.size(); ++kmer_i) {
      // Protect access to kmerMap_ts with a mutex
      //std::lock_guard<std::mutex> lock(mtx_t);
      //fill kmerMap_t with kmers and link them to sequence t id
      //kmerMap_t[seq_t_i_kmers_sorted[kmer_i]].push_back(t_i);
      //kmerMap_t_i[seq_t_i_kmers_sorted[kmer_i]]=t_i;
    //}
    kmerMap_ts[t_i] = kmerMap_t_i;
  }, ncores);
  //combine maps
  kmerMap_t = combineInt64Maps(kmerMap_ts);
  if (debug) {
    std::cout << "kmerMap_q size before: " << kmerMap_q.size();
    std::cout << std::endl;
    std::cout << "kmerMap_t size before: " << kmerMap_t.size();
    std::cout << std::endl;
  }
  //direct jaccard
  std::vector<std::vector<double>> jaccard_distances;
  jaccard_distances = directJaccardInt64(
    kmerMap_q,
    kmerMap_t,
    comparisonResults,
    seq_q_kmers_sorted,
    seq_t_kmers_sorted,
    seq_q_kmers_counts_sorted,
    seq_t_kmers_counts_sorted,
    k,
    min_jaccard);
  /*
  //create vector to store common_kmers and remove from the maps
  std::vector<std::string> common_kmers_q_t;
  common_kmers_q_t = removeUncommonKeys(kmerMap_q, kmerMap_t);
  if (debug) {
    std::cout << "kmerMap_q size after: " << kmerMap_q.size();
    std::cout << std::endl;
    std::cout << "kmerMap_t size after: " << kmerMap_t.size();
    std::cout << std::endl;
  }
  //std::set<std::pair<int, int>> pairs_set;
  //pairs_set = getPairs(kmerMap_q, kmerMap_t);
  //std::vector<std::pair<int, int>> pairs(pairs_set.begin(), pairs_set.end());
  //removeDuplicatePairs(pairs);
  std::mutex mtx_pairs;
  std::vector<std::pair<int, int>> pairs;
  RcppThread::parallelFor(0, common_kmers_q_t.size(), [&] (int kmer_i) {
    const std::string& key = common_kmers_q_t[kmer_i];
    for (const auto& pq : kmerMap_q[key]) {
      for (const auto& pt : kmerMap_t[key]) {
        std::lock_guard<std::mutex> lock(mtx_pairs);
        pairs.push_back(std::make_pair(pq, pt));
      }
    }
  }, ncores);
  if (debug) {
    std::cout << "pairs size before duplicate remove: " << pairs.size();
    std::cout << std::endl;
  }
  removeDuplicatePairs(pairs);
  if (debug) {
    std::cout << "pairs size after duplicate remove: " << pairs.size();
    std::cout << std::endl;
  }
  std::vector<std::vector<double>> jaccard_distances(pairs.size());
  jaccard_distances.resize(pairs.size());
  RcppThread::parallelFor(0, pairs.size(), [&] (int pair_i) {
    jaccard_distances[pair_i] =
      getJaccardByPair(
        pairs[pair_i],
        seq_q_kmers_sorted,
        seq_t_kmers_sorted,
        seq_q_kmers_counts_sorted,
        seq_t_kmers_counts_sorted,
        k);
  }, ncores);
  //print jaccard distances to console
  if (debug) {
    for (const auto& row : jaccard_distances) {
      for (const auto& element : row) {
        std::cout << element << " ";
      }
      std::cout << std::endl;
    }
  }
   */
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
