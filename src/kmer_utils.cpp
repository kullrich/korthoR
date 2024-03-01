//#include <Rcpp.h>
#include <RcppArmadillo.h>
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
#include "kmer_utils.h"

bool operator==(const std::tuple<std::string, std::string>& lhs, const std::tuple<std::string, std::string>& rhs) {
  return std::tie(std::get<0>(lhs), std::get<1>(lhs)) == std::tie(std::get<0>(rhs), std::get<1>(rhs));
}

std::tuple<std::vector<std::string>, std::vector<int>> sort_string_vec_and_int_vec(
  const std::vector<std::string>& string_vec,
  const std::vector<int>& int_vec) {
  std::vector<int> indices(string_vec.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
    return string_vec[i] < string_vec[j];
  });
  std::vector<std::string> sorted_string_vec(string_vec.size());
  for (int i = 0; i < string_vec.size(); ++i) {
    sorted_string_vec[i] = string_vec[indices[i]];
  }
  std::vector<int> sorted_int_vec(int_vec.size());
  for (int i = 0; i < int_vec.size(); ++i) {
    sorted_int_vec[i] = int_vec[indices[i]];
  }
  return std::make_tuple(sorted_string_vec, sorted_int_vec);
}

std::tuple<std::vector<std::int64_t>, std::vector<int>> sort_string_vec_and_int_vecInt64(
  const std::vector<std::int64_t>& string_vec,
  const std::vector<int>& int_vec) {
  std::vector<int> indices(string_vec.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int i, int j) {
    return string_vec[i] < string_vec[j];
  });
  std::vector<std::int64_t> sorted_string_vec(string_vec.size());
  for (int i = 0; i < string_vec.size(); ++i) {
    sorted_string_vec[i] = string_vec[indices[i]];
  }
  std::vector<int> sorted_int_vec(int_vec.size());
  for (int i = 0; i < int_vec.size(); ++i) {
    sorted_int_vec[i] = int_vec[indices[i]];
  }
  return std::make_tuple(sorted_string_vec, sorted_int_vec);
}

std::vector<std::string> convertRcppStringVectorToStdVector(
  const Rcpp::StringVector& string_vec) {
  std::vector<std::string> result;
  for (int i = 0; i < string_vec.size(); ++i) {
    result.push_back(Rcpp::as<std::string>(string_vec[i]));
  }
  return result;
}

std::vector<int> convertRcppIntegerVectorToStdVector(
  const Rcpp::IntegerVector& int_vec) {
  std::vector<int> result;
  for (int i = 0; i < int_vec.size(); ++i) {
    int int_vec_i;
    int_vec_i = int_vec[i];
    result.push_back(int_vec_i);
  }
  return result;
}

std::vector<std::vector<int>> convertRcppListOfIntegerVectorToStdVector(
  Rcpp::List l) {
  int n = l.size();
  std::vector<std::vector<int>> result(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::IntegerVector vec = l[i];
    result[i] = Rcpp::as<std::vector<int>>(vec);
  }
  return result;
}

/*
std::vector<std::vector<std::string>> convertInnerNamesRcppListOfIntegerVectorToStdVector(
  Rcpp::List l) {
  int n = l.size();
  std::vector<std::vector<std::string>> result(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::IntegerVector vec = l[i];
    std::vector<std::string> vecnames = vec.attr("names");
    result[i] = vecnames;
  }
  return result;
}
*/

std::vector<std::vector<std::string>> convertInnerNamesRcppListOfIntegerVectorToStdVector(
  Rcpp::List l) {
  int n = l.size();
  std::vector<std::vector<std::string>> result(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::RObject obj = l[i];
    SEXP objNames = obj.attr("names");
    if(!Rf_isNull(objNames)) {
      Rcpp::CharacterVector index(objNames);
      std::vector<std::string> names(index.size());
      for(int j = 0; j < index.size(); ++j) {
        names[j] = Rcpp::as<std::string>(index[j]);
      }
      result[i] = names;
    }
  }
  return result;
}

std::vector<std::vector<std::int64_t>> convertInnerNamesRcppListOfIntegerVectorToStdVectorInt64(
  Rcpp::List l) {
  int n = l.size();
  std::vector<std::vector<std::int64_t>> result(n);
  for (int i = 0; i < n; ++i) {
    Rcpp::RObject obj = l[i];
    SEXP objNames = obj.attr("names");
    if(!Rf_isNull(objNames)) {
      Rcpp::CharacterVector index(objNames);
      std::vector<std::int64_t> names(index.size());
      for(int j = 0; j < index.size(); ++j) {
        std::string index_j;
        index_j = index[j];
        names[j] = std::stoll(index_j);
      }
      result[i] = names;
    }
  }
  return result;
}

Rcpp::IntegerVector name_int_vec_by_string_vec(
  const std::vector<int>& int_vec,
  const std::vector<std::string>& string_vec) {
  Rcpp::IntegerVector result_int_vec(int_vec.size());
  Rcpp::CharacterVector result_char_vec(string_vec.size());
  for (int i = 0; i < int_vec.size(); ++i) {
    result_int_vec[i] = int_vec[i];
  }
  for (int j = 0; j < string_vec.size(); ++j) {
    result_char_vec[j] = string_vec[j];
  }
  result_int_vec.attr("names") = result_char_vec;
  return result_int_vec;
}

Rcpp::IntegerVector name_int_vec_by_string_vecInt64(
  const std::vector<int>& int_vec,
  const std::vector<std::int64_t>& string_vec) {
  Rcpp::IntegerVector result_int_vec(int_vec.size());
  Rcpp::NumericVector result_char_vec(string_vec.size());
  for (int i = 0; i < int_vec.size(); ++i) {
    result_int_vec[i] = int_vec[i];
  }
  for (int j = 0; j < string_vec.size(); ++j) {
    result_char_vec[j] = string_vec[j];
  }
  result_int_vec.attr("names") = result_char_vec;
  return result_int_vec;
}

Rcpp::List iterate_and_create_named_int_vec(
  const std::vector<std::vector<int>> int_vec_array,
  const std::vector<std::vector<std::string>> string_vec_array,
  const std::vector<std::string> seqnames) {
  Rcpp::List result_int_vec(int_vec_array.size());
  for (int i = 0; i < int_vec_array.size(); ++i) {
    result_int_vec[i] = name_int_vec_by_string_vec(int_vec_array[i], string_vec_array[i]);
  }
  result_int_vec.attr("names") = seqnames;
  return result_int_vec;
}

Rcpp::List iterate_and_create_named_int_vecInt64(
  const std::vector<std::vector<int>> int_vec_array,
  const std::vector<std::vector<std::int64_t>> string_vec_array,
  const std::vector<std::string> seqnames) {
  Rcpp::List result_int_vec(int_vec_array.size());
  for (int i = 0; i < int_vec_array.size(); ++i) {
    result_int_vec[i] = name_int_vec_by_string_vecInt64(int_vec_array[i], string_vec_array[i]);
  }
  result_int_vec.attr("names") = seqnames;
  return result_int_vec;
}

std::map<std::string, int> vectorToMap(
  const std::vector<std::string>& vec,
  int value) {
  std::map<std::string, int> resultMap;
  for (const std::string& str : vec) {
    resultMap[str] = value;
  }
  return resultMap;
}

std::map<std::int64_t, int> vectorToInt64Map(
  const std::vector<std::int64_t>& vec,
  int value) {
  std::map<std::int64_t, int> resultInt64Map;
  for (const std::int64_t& str : vec) {
    resultInt64Map[str] = value;
  }
  return resultInt64Map;
}

std::map<std::string, std::vector<int>> combineMaps(
  const std::vector<std::map<std::string, int>>& maps) {
    std::map<std::string, std::vector<int>> combinedMap;
    for (const auto& map : maps) {
      for (const auto& pair : map) {
        combinedMap[pair.first].push_back(pair.second);
    }
  }
  return combinedMap;
}

std::map<std::string, std::vector<std::vector<int>>> combineMapsTwo(
  const std::vector<std::map<std::string, int>>& maps1,
  const std::vector<std::map<std::string, int>>& maps2) {
  std::map<std::string, std::vector<std::vector<int>>> combinedMap;
  for (const std::map<std::string, int>& map1 : maps1) {
    for (const auto& pair : map1) {
      combinedMap[pair.first].resize(2);
      for (const std::map<std::string, int>& map2 : maps2) {
        try {
          combinedMap[pair.first][1].push_back(map2.at(pair.first));
        } catch (const std::out_of_range& e) {
        }
      }
      if (combinedMap[pair.first][1].size() > 0) {
        combinedMap[pair.first][0].push_back(pair.second);
      } else {
        combinedMap.erase(pair.first);
      }
    }
  }
  return combinedMap;
}

std::map<std::int64_t, std::vector<int>> combineInt64Maps(
  const std::vector<std::map<std::int64_t, int>>& maps) {
    std::map<std::int64_t, std::vector<int>> combinedInt64Map;
    for (const auto& map : maps) {
      for (const auto& pair : map) {
        combinedInt64Map[pair.first].push_back(pair.second);
    }
  }
  return combinedInt64Map;
}

std::vector<int> CountCommonAndUncommonAndSum(
  const std::vector<std::string>& vec1,
  const std::vector<std::string>& vec2,
  const std::vector<int>& vec1_counts,
  const std::vector<int>& vec2_counts) {
  int i = 0;
  int j = 0;
  int common = 0;
  int uncommon1 = 0;
  int uncommon2 = 0;
  int common1_sum = 0;
  int common2_sum = 0;
  int uncommon1_sum = 0;
  int uncommon2_sum = 0;
  std::vector<int> result(7);
  while (i < vec1.size() && j < vec2.size()) {
    if (vec1[i] < vec2[j]) {
      uncommon1 += 1;
      uncommon1_sum += vec1_counts[i];
      ++i;
    } else if (vec1[i] > vec2[j]) {
      uncommon2 += 1;
      uncommon2_sum += vec2_counts[j];
      ++j;
    } else {
      common += 1;
      common1_sum += vec1_counts[i];
      common2_sum += vec2_counts[j];
      ++i;
      ++j;
    }
  }
  for (; i < vec1.size(); ++i) {
    uncommon1 += 1;
    uncommon1_sum += vec1_counts[i];
  }
  for (; j < vec2.size(); ++j) {
    uncommon2 += 1;
    uncommon2_sum += vec2_counts[j];
  }
  result[0] = common;
  result[1] = uncommon1;
  result[2] = uncommon2;
  result[3] = common1_sum;
  result[4] = common2_sum;
  result[5] = uncommon1_sum;
  result[6] = uncommon2_sum;
  return result;
}

std::vector<double> get_jaccard(
  const int common,
  const int uncommon1,
  const int uncommon2,
  const int k) {
  std::vector<double> seqj_seqk_jaccard(3);
  double jaccard;
  double mash;
  double ani;
  jaccard = ( static_cast<double>(common) ) / ( common + uncommon1 + uncommon2 );
  mash = -( 1/static_cast<double>(k) ) * std::log( ( 2 * jaccard ) / ( 1 + jaccard ) );
  ani = 1 - mash;
  seqj_seqk_jaccard[0] = jaccard;
  seqj_seqk_jaccard[1] = mash;
  seqj_seqk_jaccard[2] = ani;
  return seqj_seqk_jaccard;
}

double get_jaccard_value(
  const int common,
  const int uncommon1,
  const int uncommon2) {
  double jaccard_value;
  jaccard_value = ( static_cast<double>(common) ) / ( common + uncommon1 + uncommon2 );
  return jaccard_value;
}

double get_mash_value(
  const double jaccard_value,
  const int k) {
  double mash_value;
  mash_value = -( 1/static_cast<double>(k) ) * std::log( ( 2 * jaccard_value ) / ( 1 + jaccard_value ) );
  return mash_value;
}

double get_ani_value(
  const double mash_value) {
  double ani_value;
  ani_value = 1 - mash_value;
  return ani_value;
}

double get_sumdist_value(
  const int common1_sum,
  const int common2_sum,
  const int uncommon1_sum,
  const int uncommon2_sum) {
  double sumdist_value;
  sumdist_value = 1-( static_cast<double>(common1_sum) + common2_sum ) / ( common1_sum + common2_sum + uncommon1_sum + uncommon2_sum );
  return sumdist_value;
}

std::vector<std::vector<double>> directJaccard(
  std::map<std::string, std::vector<int>>& map1,
  std::map<std::string, std::vector<int>>& map2,
  std::vector<std::vector<bool>> comparisonResults,
  const std::vector<std::vector<std::string>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::string>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k,
  const double min_jaccard,
  const int ncores) {
  //std::mutex mtx_distances;
  std::vector<std::vector<double>> jaccard_distances;
  auto it1 = map1.begin();
  auto it2 = map2.begin();
  while (it1 != map1.end() && it2 != map2.end()) {
    if (it1->first < it2->first) {
      it1++;
    } else if (it2->first < it1->first) {
      it2++;
    } else {
      if (it1->second.size() > it2->second.size()) {
        for (size_t q_i = 0; q_i < it1->second.size(); ++q_i) {
          for (size_t t_j = 0; t_j < it2->second.size(); ++t_j) {
            if (comparisonResults[it1->second[q_i]][it2->second[t_j]]) {
            } else {
              std::vector<double> jaccard_distance_q_i_t_j(6);
              jaccard_distance_q_i_t_j = getJaccardByPair(
                std::make_pair(it1->second[q_i], it2->second[t_j]),
                seq_q_kmers_sorted,
                seq_t_kmers_sorted,
                seq_q_kmers_counts_sorted,
                seq_t_kmers_counts_sorted,
                k);
              if (jaccard_distance_q_i_t_j[2] > min_jaccard) {
                //std::lock_guard<std::mutex> lock(mtx_distances);
                jaccard_distances.push_back(jaccard_distance_q_i_t_j);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              } else {
                //std::lock_guard<std::mutex> lock(mtx_distances);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              }
            }
          }
        }
      } else {
        for (size_t t_j = 0; t_j < it2->second.size(); ++t_j) {
          for (size_t q_i = 0; q_i < it1->second.size(); ++q_i) {
            if (comparisonResults[it1->second[q_i]][it2->second[t_j]]) {
            } else {
              std::vector<double> jaccard_distance_q_i_t_j(6);
              jaccard_distance_q_i_t_j = getJaccardByPair(
                std::make_pair(it1->second[q_i], it2->second[t_j]),
                seq_q_kmers_sorted,
                seq_t_kmers_sorted,
                seq_q_kmers_counts_sorted,
                seq_t_kmers_counts_sorted,
                k);
              if (jaccard_distance_q_i_t_j[2] > min_jaccard) {
                //std::lock_guard<std::mutex> lock(mtx_distances);
                jaccard_distances.push_back(jaccard_distance_q_i_t_j);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              } else {
                //std::lock_guard<std::mutex> lock(mtx_distances);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              }
            }
          }
        }
      }
      ++it1;
      ++it2;
    }
  }
  // Not need to process further
  return jaccard_distances;
}

std::vector<std::vector<double>> directJaccardInt64(
  std::map<std::int64_t, std::vector<int>>& map1,
  std::map<std::int64_t, std::vector<int>>& map2,
  std::vector<std::vector<bool>> comparisonResults,
  const std::vector<std::vector<std::int64_t>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::int64_t>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k,
  const double min_jaccard,
  const int ncores) {
  std::mutex mtx_distances;
  std::vector<std::vector<double>> jaccard_distances;
  auto it1 = map1.begin();
  auto it2 = map2.begin();
  while (it1 != map1.end() && it2 != map2.end()) {
    if (it1->first < it2->first) {
      it1++;
    } else if (it2->first < it1->first) {
      it2++;
    } else {
      if (it1->second.size() > it2->second.size()) {
        for (size_t q_i = 0; q_i < it1->second.size(); ++q_i) {
          for (size_t t_j = 0; t_j < it2->second.size(); ++t_j) {
            if (comparisonResults[it1->second[q_i]][it2->second[t_j]]) {
            } else {
              std::vector<double> jaccard_distance_q_i_t_j(6);
              jaccard_distance_q_i_t_j = getJaccardByPairInt64(
                std::make_pair(it1->second[q_i], it2->second[t_j]),
                seq_q_kmers_sorted,
                seq_t_kmers_sorted,
                seq_q_kmers_counts_sorted,
                seq_t_kmers_counts_sorted,
                k);
              if (jaccard_distance_q_i_t_j[2] > min_jaccard) {
                std::lock_guard<std::mutex> lock(mtx_distances);
                jaccard_distances.push_back(jaccard_distance_q_i_t_j);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              } else {
                std::lock_guard<std::mutex> lock(mtx_distances);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              }
            }
          }
        }
      } else {
        for (size_t t_j = 0; t_j < it2->second.size(); ++t_j) {
          for (size_t q_i = 0; q_i < it1->second.size(); ++q_i) {
            if (comparisonResults[it1->second[q_i]][it2->second[t_j]]) {
            } else {
              std::vector<double> jaccard_distance_q_i_t_j(6);
              jaccard_distance_q_i_t_j = getJaccardByPairInt64(
                std::make_pair(it1->second[q_i], it2->second[t_j]),
                seq_q_kmers_sorted,
                seq_t_kmers_sorted,
                seq_q_kmers_counts_sorted,
                seq_t_kmers_counts_sorted,
                k);
              if (jaccard_distance_q_i_t_j[2] > min_jaccard) {
                std::lock_guard<std::mutex> lock(mtx_distances);
                jaccard_distances.push_back(jaccard_distance_q_i_t_j);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              } else {
                std::lock_guard<std::mutex> lock(mtx_distances);
                comparisonResults[it1->second[q_i]][it2->second[t_j]] = true;
              }
            }
          }
        }
      }
      ++it1;
      ++it2;
    }
  }
  // Not need to process further
  return jaccard_distances;
}

std::vector<std::string> removeUncommonKeys(
  std::map<std::string,
  std::vector<int>>& map1,
  std::map<std::string,
  std::vector<int>>& map2) {
  std::vector<std::string> common_kmers_q_t;
  common_kmers_q_t.reserve(std::min(map1.size(), map2.size()));
  std::vector<std::string> keys_to_erase_map1;
  std::vector<std::string> keys_to_erase_map2;
  auto it1 = map1.begin();
  auto it2 = map2.begin();
  while (it1 != map1.end() && it2 != map2.end()) {
    if (it1->first < it2->first) {
      keys_to_erase_map1.push_back(it1->first);
      ++it1;
    } else if (it2->first < it1->first) {
      keys_to_erase_map2.push_back(it2->first);
      ++it2;
    } else {
      common_kmers_q_t.push_back(it1->first);
      ++it1;
      ++it2;
    }
  }
  // Remove remaining keys in map1
  for (; it1 != map1.end(); ++it1) {
    keys_to_erase_map1.push_back(it1->first);
  }
  // Remove remaining keys in map2
  for (; it2 != map2.end(); ++it2) {
    keys_to_erase_map2.push_back(it2->first);
  }
  // Erase marked keys
  for (const auto& key : keys_to_erase_map1) {
    map1.erase(key);
  }
  for (const auto& key : keys_to_erase_map2) {
    map2.erase(key);
  }
  return common_kmers_q_t;
}

std::set<std::pair<int, int>> getPairs(
  std::map<std::string, std::vector<int>>& map1,
  std::map<std::string, std::vector<int>>& map2) {
  std::set<std::pair<int, int>> pairs;
  auto it1 = map1.begin();
  auto it2 = map2.begin();
  while (it1 != map1.end() && it2 != map2.end()) {
    if (it1->first < it2->first) {
    } else if (it2->first < it1->first) {
    } else {
      for (const auto& pq : it1->second) {
        for (const auto& pt : it2->second) {
          pairs.insert(std::make_pair(pq, pt));
        }
      }
      ++it1;
      ++it2;
    }
  }
  return pairs;
}

void removeDuplicatePairs(
    std::vector<std::pair<int, int>>& pairs) {
    std::sort(pairs.begin(), pairs.end());
    auto last = std::unique(pairs.begin(), pairs.end());
    pairs.erase(last, pairs.end());
}

std::vector<double> getJaccardByIntegerVector(
  const std::vector<std::string> seq_q_i_kmers_sorted,
  const std::vector<std::string> seq_t_i_kmers_sorted,
  const std::vector<int> seq_q_i_kmers_counts_sorted,
  const std::vector<int> seq_t_i_kmers_counts_sorted,
  const int seq_q_i_idx,
  const int seq_t_i_idx,
  const int k) {
  std::vector<double> out(6);
  out[0] = static_cast<double>(seq_q_i_idx);
  out[1] = static_cast<double>(seq_t_i_idx);
  int i = 0;
  int j = 0;
  int common = 0;
  int uncommon1 = 0;
  int uncommon2 = 0;
  int common1_sum = 0;
  int common2_sum = 0;
  int uncommon1_sum = 0;
  int uncommon2_sum = 0;
  while (i < seq_q_i_kmers_sorted.size() && j < seq_t_i_kmers_sorted.size()) {
    if (seq_q_i_kmers_sorted[i] < seq_t_i_kmers_sorted[j]) {
      uncommon1 += 1;
      uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
      ++i;
    } else if (seq_q_i_kmers_sorted[i] > seq_t_i_kmers_sorted[j]) {
      uncommon2 += 1;
      uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
      ++j;
    } else {
      common += 1;
      common1_sum += seq_q_i_kmers_counts_sorted[i];
      common2_sum += seq_t_i_kmers_counts_sorted[j];
      ++i;
      ++j;
    }
  }
  for (; i < seq_q_i_kmers_sorted.size(); ++i) {
    uncommon1 += 1;
    uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
  }
  for (; j < seq_t_i_kmers_sorted.size(); ++j) {
    uncommon2 += 1;
    uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
  }
  //calculate jaccard, mash, ani, sumdist
  double jaccard_value;
  double mash_value;
  double ani_value;
  double sumdist_value;
  jaccard_value = get_jaccard_value(common, uncommon1, uncommon2);
  mash_value = get_mash_value(jaccard_value, k);
  ani_value = get_ani_value(mash_value);
  sumdist_value = get_sumdist_value(common1_sum, common2_sum, uncommon1_sum, uncommon2_sum);
  out[2] = jaccard_value;
  out[3] = mash_value;
  out[4] = ani_value;
  out[5] = sumdist_value;
  return out;
}

std::vector<double> getJaccardByIntegerVectorInt64(
  const std::vector<std::int64_t> seq_q_i_kmers_sorted,
  const std::vector<std::int64_t> seq_t_i_kmers_sorted,
  const std::vector<int> seq_q_i_kmers_counts_sorted,
  const std::vector<int> seq_t_i_kmers_counts_sorted,
  const int seq_q_i_idx,
  const int seq_t_i_idx,
  const int k) {
  std::vector<double> out(6);
  out[0] = static_cast<double>(seq_q_i_idx);
  out[1] = static_cast<double>(seq_t_i_idx);
  int i = 0;
  int j = 0;
  int common = 0;
  int uncommon1 = 0;
  int uncommon2 = 0;
  int common1_sum = 0;
  int common2_sum = 0;
  int uncommon1_sum = 0;
  int uncommon2_sum = 0;
  while (i < seq_q_i_kmers_sorted.size() && j < seq_t_i_kmers_sorted.size()) {
    if (seq_q_i_kmers_sorted[i] < seq_t_i_kmers_sorted[j]) {
      uncommon1 += 1;
      uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
      ++i;
    } else if (seq_q_i_kmers_sorted[i] > seq_t_i_kmers_sorted[j]) {
      uncommon2 += 1;
      uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
      ++j;
    } else {
      common += 1;
      common1_sum += seq_q_i_kmers_counts_sorted[i];
      common2_sum += seq_t_i_kmers_counts_sorted[j];
      ++i;
      ++j;
    }
  }
  for (; i < seq_q_i_kmers_sorted.size(); ++i) {
    uncommon1 += 1;
    uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
  }
  for (; j < seq_t_i_kmers_sorted.size(); ++j) {
    uncommon2 += 1;
    uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
  }
  //calculate jaccard, mash, ani, sumdist
  double jaccard_value;
  double mash_value;
  double ani_value;
  double sumdist_value;
  jaccard_value = get_jaccard_value(common, uncommon1, uncommon2);
  mash_value = get_mash_value(jaccard_value, k);
  ani_value = get_ani_value(mash_value);
  sumdist_value = get_sumdist_value(common1_sum, common2_sum, uncommon1_sum, uncommon2_sum);
  out[2] = jaccard_value;
  out[3] = mash_value;
  out[4] = ani_value;
  out[5] = sumdist_value;
  return out;
}

std::vector<double> getJaccardByPair(
  const std::pair<int, int>& p,
  const std::vector<std::vector<std::string>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::string>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k) {
  std::vector<double> out(6);
  out[0] = static_cast<double>(p.first);
  out[1] = static_cast<double>(p.second);
  int i = 0;
  int j = 0;
  int common = 0;
  int uncommon1 = 0;
  int uncommon2 = 0;
  int common1_sum = 0;
  int common2_sum = 0;
  int uncommon1_sum = 0;
  int uncommon2_sum = 0;
  //get pair kmers and kmers counts
  std::vector<std::string> seq_q_i_kmers_sorted = seq_q_kmers_sorted[p.first];
  std::vector<std::string> seq_t_i_kmers_sorted = seq_t_kmers_sorted[p.second];
  std::vector<int> seq_q_i_kmers_counts_sorted = seq_q_kmers_counts_sorted[p.first];
  std::vector<int> seq_t_i_kmers_counts_sorted = seq_t_kmers_counts_sorted[p.second];
  while (i < seq_q_i_kmers_sorted.size() && j < seq_t_i_kmers_sorted.size()) {
    if (seq_q_i_kmers_sorted[i] < seq_t_i_kmers_sorted[j]) {
      uncommon1 += 1;
      uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
      ++i;
    } else if (seq_q_i_kmers_sorted[i] > seq_t_i_kmers_sorted[j]) {
      uncommon2 += 1;
      uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
      ++j;
    } else {
      common += 1;
      common1_sum += seq_q_i_kmers_counts_sorted[i];
      common2_sum += seq_t_i_kmers_counts_sorted[j];
      ++i;
      ++j;
    }
  }
  for (; i < seq_q_i_kmers_sorted.size(); ++i) {
    uncommon1 += 1;
    uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
  }
  for (; j < seq_t_i_kmers_sorted.size(); ++j) {
    uncommon2 += 1;
    uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
  }
  //calculate jaccard, mash, ani, sumdist
  double jaccard_value;
  double mash_value;
  double ani_value;
  double sumdist_value;
  jaccard_value = get_jaccard_value(common, uncommon1, uncommon2);
  mash_value = get_mash_value(jaccard_value, k);
  ani_value = get_ani_value(mash_value);
  sumdist_value = get_sumdist_value(common1_sum, common2_sum, uncommon1_sum, uncommon2_sum);
  out[2] = jaccard_value;
  out[3] = mash_value;
  out[4] = ani_value;
  out[5] = sumdist_value;
  return out;
}

std::vector<double> getJaccardByPairInt64(
  const std::pair<int, int>& p,
  const std::vector<std::vector<std::int64_t>>& seq_q_kmers_sorted,
  const std::vector<std::vector<std::int64_t>>& seq_t_kmers_sorted,
  const std::vector<std::vector<int>>& seq_q_kmers_counts_sorted,
  const std::vector<std::vector<int>>& seq_t_kmers_counts_sorted,
  const int k) {
  std::vector<double> out(6);
  out[0] = static_cast<double>(p.first);
  out[1] = static_cast<double>(p.second);
  int i = 0;
  int j = 0;
  int common = 0;
  int uncommon1 = 0;
  int uncommon2 = 0;
  int common1_sum = 0;
  int common2_sum = 0;
  int uncommon1_sum = 0;
  int uncommon2_sum = 0;
  //get pair kmers and kmers counts
  std::vector<std::int64_t> seq_q_i_kmers_sorted = seq_q_kmers_sorted[p.first];
  std::vector<std::int64_t> seq_t_i_kmers_sorted = seq_t_kmers_sorted[p.second];
  std::vector<int> seq_q_i_kmers_counts_sorted = seq_q_kmers_counts_sorted[p.first];
  std::vector<int> seq_t_i_kmers_counts_sorted = seq_t_kmers_counts_sorted[p.second];
  while (i < seq_q_i_kmers_sorted.size() && j < seq_t_i_kmers_sorted.size()) {
    if (seq_q_i_kmers_sorted[i] < seq_t_i_kmers_sorted[j]) {
      uncommon1 += 1;
      uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
      ++i;
    } else if (seq_q_i_kmers_sorted[i] > seq_t_i_kmers_sorted[j]) {
      uncommon2 += 1;
      uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
      ++j;
    } else {
      common += 1;
      common1_sum += seq_q_i_kmers_counts_sorted[i];
      common2_sum += seq_t_i_kmers_counts_sorted[j];
      ++i;
      ++j;
    }
  }
  for (; i < seq_q_i_kmers_sorted.size(); ++i) {
    uncommon1 += 1;
    uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
  }
  for (; j < seq_t_i_kmers_sorted.size(); ++j) {
    uncommon2 += 1;
    uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
  }
  //calculate jaccard, mash, ani, sumdist
  double jaccard_value;
  double mash_value;
  double ani_value;
  double sumdist_value;
  jaccard_value = get_jaccard_value(common, uncommon1, uncommon2);
  mash_value = get_mash_value(jaccard_value, k);
  ani_value = get_ani_value(mash_value);
  sumdist_value = get_sumdist_value(common1_sum, common2_sum, uncommon1_sum, uncommon2_sum);
  out[2] = jaccard_value;
  out[3] = mash_value;
  out[4] = ani_value;
  out[5] = sumdist_value;
  return out;
}

std::int64_t aminoAcidSequenceToInt(
  const std::string& sequence) {
  std::string prefixedSequence = "#" + sequence;
  auto aminoAcidToInt = [](char aa) -> int {
    switch (aa) {
      case 'A': return 0;
      case 'C': return 1;
      case 'D': return 2;
      case 'E': return 3;
      case 'F': return 4;
      case 'G': return 5;
      case 'H': return 6;
      case 'I': return 7;
      case 'K': return 8;
      case 'L': return 9;
      case 'M': return 10;
      case 'N': return 11;
      case 'P': return 12;
      case 'Q': return 13;
      case 'R': return 14;
      case 'S': return 15;
      case 'T': return 16;
      case 'V': return 17;
      case 'W': return 18;
      case 'Y': return 19;
      case '*': return 20;
      default: return 21; // Handle unknown amino acids
    }
  };
  std::int64_t numericRepresentation = 0;
  for (char aa : prefixedSequence) {
    numericRepresentation = numericRepresentation * 22 + aminoAcidToInt(aa);
  }
  return numericRepresentation;
}

arma::imat changeValueAndUpdateSparseImat(
  const arma::imat& sp_mat_logical,
  int row,
  int col,
  int newValue) {
  arma::imat modified_mat = sp_mat_logical;
  modified_mat(row, col) = newValue;
  return modified_mat;
}

bool checkValueSparseImat(
  const arma::imat& sp_mat_logical,
  int row,
  int col) {
  return sp_mat_logical(row, col) == 1;
}

bool key_found(
  const std::map<std::string, int>& small_map,
  const std::map<std::string, int>& large_map) {
  for (const auto& entry : small_map) {
    auto it = large_map.lower_bound(entry.first);
    if (it != large_map.end() && !(large_map.key_comp()(entry.first, it->first))) {
      return true;
    }
  }
  return false;
}

bool key_found2(
  const std::map<std::string, int>& small_map,
  const std::map<std::string, int>& large_map) {
  auto it1 = small_map.begin();
  auto it2 = large_map.begin();
  while (it1 != small_map.end() && it2 != large_map.end()) {
    if (it1->first < it2->first) {
      ++it1;
    } else if (it2->first < it1->first) {
      ++it2;
    } else {
      return true;
      ++it1;
      ++it2;
    }
  }
  return false;
}

std::vector<int> findHits(
  const std::map<std::string, int>& kmerMap_i,
  const std::map<std::string, std::vector<int>>& kmerMap) {
  auto it1 = kmerMap_i.begin();
  auto it2 = kmerMap.begin();
  std::vector<int> kmerMap_i_hits;
  while (it1 != kmerMap_i.end() && it2 != kmerMap.end()) {
    if (it1->first < it2->first) {
      ++it1;
    } else if (it2->first < it1->first) {
      ++it2;
    } else {
      for (const int& hit : it2->second) {
        kmerMap_i_hits.push_back(hit);
      }
      ++it1;
      ++it2;
    }
  }
  std::set<int> kmerMap_i_hits_set(kmerMap_i_hits.begin(), kmerMap_i_hits.end());
  std::vector<int> result(kmerMap_i_hits_set.size());
  result = std::vector<int>(kmerMap_i_hits_set.begin(), kmerMap_i_hits_set.end());
  return result;
}

std::vector<std::vector<int>> findHits_sparse(
  const std::map<std::string, std::vector<int>>& kmerMap_short,
  const std::map<std::string, std::vector<int>>& kmerMap_long,
  const int n) {
  auto it1 = kmerMap_short.begin();
  auto it2 = kmerMap_long.begin();
  std::vector<std::vector<int>> kmerMap_i_j_hits(n);
  while (it1 != kmerMap_short.end() && it2 != kmerMap_long.end()) {
    if (it1->first < it2->first) {
      ++it1;
    } else if (it2->first < it1->first) {
      ++it2;
    } else {
      for (const int& i : it1->second) {
        for (const int& j : it2->second) {
          kmerMap_i_j_hits[i].push_back(j);
        }
      }
      ++it1;
      ++it2;
    }
  }
  std::vector<std::vector<int>> result(n);
  for (int i = 0; i < n; ++i) {
    std::set<int> kmerMap_i_j_hits_set(kmerMap_i_j_hits[i].begin(), kmerMap_i_j_hits[i].end());
    std::vector<int> result_i(kmerMap_i_j_hits_set.size());
    result[i] = std::vector<int>(kmerMap_i_j_hits_set.begin(), kmerMap_i_j_hits_set.end());
  }
  return result;
}

std::map<std::string, std::vector<int>> getSubset(
    const std::map<std::string, std::vector<int>>& originalMap,
    int n) {
    std::map<std::string, std::vector<int>> subsetMap;
    int count = 0;
    for (const auto& pair : originalMap) {
        if (count >= n) {
            break;
        }
        subsetMap.insert(pair);
        count++;
    }
    return subsetMap;
}

std::vector<std::pair<int, int>> flatten(
  const std::vector<std::vector<int>>& vec) {
  std::vector<std::pair<int, int>> result;
  for (int i = 0; i < vec.size(); ++i) {
    for (int j = 0; j < vec[i].size(); ++j) {
      result.push_back(std::make_pair(i, vec[i][j]));
    }
  }
  return result;
}
