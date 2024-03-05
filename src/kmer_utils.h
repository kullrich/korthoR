#ifndef kmer_utils_h
#define kmer_utils_h
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override

#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iterator>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace std {
  template <>
  struct hash<std::tuple<std::string, std::string>> {
    size_t operator()(const std::tuple<std::string, std::string>& tuple) const {
      size_t hashValue = 0;
      hashValue ^= std::hash<std::string>{}(std::get<0>(tuple));
      hashValue ^= std::hash<std::string>{}(std::get<1>(tuple)) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
      return hashValue;
    }
  };
}

std::tuple<std::vector<std::string>, std::vector<int>> sort_string_vec_and_int_vec(
  const std::vector<std::string>& string_vec,
  const std::vector<int>& int_vec);
std::tuple<std::vector<std::int64_t>, std::vector<int>> sort_string_vec_and_int_vecInt64(
  const std::vector<std::int64_t>& string_vec,
  const std::vector<int>& int_vec);
std::vector<std::string> convertRcppStringVectorToStdVector(
  const Rcpp::StringVector& string_vec);
std::vector<int> convertRcppIntegerVectorToStdVector(
  const Rcpp::IntegerVector& int_vec);
std::vector<std::vector<int>> convertRcppListOfIntegerVectorToStdVector(
  Rcpp::List l);
std::vector<std::vector<std::string>> convertInnerNamesRcppListOfIntegerVectorToStdVector(
  Rcpp::List l);
std::vector<std::vector<std::int64_t>> convertInnerNamesRcppListOfIntegerVectorToStdVectorInt64(
  Rcpp::List l);
Rcpp::IntegerVector name_int_vec_by_string_vec(
  const std::vector<int>& int_vec,
  const std::vector<std::string>& string_vec);
Rcpp::List iterate_and_create_named_int_vec(
  const std::vector<std::vector<int>> int_vec_array,
  const std::vector<std::vector<std::string>> string_vec_array,
  const std::vector<std::string> seqnames);
Rcpp::List iterate_and_create_named_int_vecInt64(
  const std::vector<std::vector<int>> int_vec_array,
  const std::vector<std::vector<std::int64_t>> string_vec_array,
  const std::vector<std::string> seqnames);
std::map<std::string, std::vector<int>> combineMaps(
  const std::vector<std::map<std::string, int>>& maps);
std::map<std::int64_t, std::vector<int>> combineInt64Maps(
  const std::vector<std::map<std::int64_t, int>>& maps);
std::map<std::string, int> vectorToMap(
  const std::vector<std::string>& vec,
  int value);
std::map<std::int64_t, int> vectorToInt64Map(
  const std::vector<std::int64_t>& vec,
  int value);
std::vector<int> CountCommonAndUncommonAndSum(
  const std::vector<std::string>& vec1,
  const std::vector<std::string>& vec2,
  const std::vector<int>& vec1_counts,
  const std::vector<int>& vec2_counts);
std::vector<double> get_jaccard(
  const int common,
  const int uncommon1,
  const int uncommon2,
  const int k);
double get_jaccard_value(
  const int common,
  const int uncommon1,
  const int uncommon2);
double get_mash_value(
  const double jaccard_value,
  const int k);
double get_ani_value(
  const double mash_value);
double get_sumdist_value(
  const int common1_sum,
  const int common2_sum,
  const int uncommon1_sum,
  const int uncommon2_sum);
std::vector<std::vector<double>> directJaccard(
  std::map<std::string, std::vector<int>>& map1,
  std::map<std::string, std::vector<int>>& map2,
  std::vector<std::vector<bool>> comparisonResults,
  const std::vector<std::vector<std::string>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::string>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k,
  const double min_jaccard);
std::vector<std::vector<double>> directJaccardInt64(
  std::map<std::int64_t, std::vector<int>>& map1,
  std::map<std::int64_t, std::vector<int>>& map2,
  std::vector<std::vector<bool>> comparisonResults,
  const std::vector<std::vector<std::int64_t>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::int64_t>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k,
  const double min_jaccard);
std::vector<std::string> removeUncommonKeys(
  std::map<std::string,
  std::vector<int>>& map1,
  std::map<std::string,
  std::vector<int>>& map2);
std::set<std::pair<int, int>> getPairs(
  std::map<std::string,
  std::vector<int>>& map1,
  std::map<std::string,
  std::vector<int>>& map2);
void removeDuplicatePairs(
    std::vector<std::pair<int, int>>& pairs);
std::vector<double> getJaccardByPair(
  const std::pair<int, int>& p,
  const std::vector<std::vector<std::string>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::string>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k);
std::vector<double> getJaccardByPairInt64(
  const std::pair<int, int>& p,
  const std::vector<std::vector<std::int64_t>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::int64_t>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k);
std::int64_t aminoAcidSequenceToInt(const std::string& sequence);
std::vector<double> getJaccardByIntegerVector(
  const std::vector<std::string> seq_q_i_kmers_sorted,
  const std::vector<std::string> seq_t_i_kmers_sorted,
  const std::vector<int> seq_q_i_kmers_counts_sorted,
  const std::vector<int> seq_t_i_kmers_counts_sorted,
  const int seq_q_i_idx,
  const int seq_t_i_idx,
  const int k);
std::vector<double> getJaccardByIntegerVectorInt64(
  const std::vector<std::int64_t> seq_q_i_kmers_sorted,
  const std::vector<std::int64_t> seq_t_i_kmers_sorted,
  const std::vector<int> seq_q_i_kmers_counts_sorted,
  const std::vector<int> seq_t_i_kmers_counts_sorted,
  const int seq_q_i_idx,
  const int seq_t_i_idx,
  const int k);
std::vector<int> findPairs(
  const std::map<std::string, int>& kmerMap_i,
  const std::map<std::string, std::vector<int>>& kmerMap);
std::vector<int> findPairsInt64(
  const std::map<std::int64_t, int>& kmerMap_i,
  const std::map<std::int64_t, std::vector<int>>& kmerMap);
std::vector<std::vector<int>> findPairs_sparse(
  const std::map<std::string, std::vector<int>>& kmerMap_short,
  const std::map<std::string, std::vector<int>>& kmerMap_long,
  const int n);
std::vector<std::vector<int>> findPairsInt64_sparse(
  const std::map<std::int64_t, std::vector<int>>& kmerMap_short,
  const std::map<std::int64_t, std::vector<int>>& kmerMap_long,
  const int n);
std::map<std::string, std::vector<int>> getSubsetMap(
    const std::map<std::string, std::vector<int>>& originalMap,
    int n);
std::map<std::int64_t, std::vector<int>> getSubsetInt64Map(
    const std::map<std::int64_t, std::vector<int>>& originalMap,
    int n);
std::vector<std::set<int>> getCandidates(
  const std::vector<std::vector<int>>& vec);
std::vector<std::pair<int, int>> flatten2pairs(
  const std::vector<std::vector<int>>& vec);
std::tuple<std::vector<int>, std::vector<int>> flatten2tuple(
  const std::vector<std::vector<int>>& vec);
std::vector<std::vector<int>> transposeTQ(
  const std::vector<std::vector<int>>& matrix,
  const int n);

#endif // kmer_utils_h
