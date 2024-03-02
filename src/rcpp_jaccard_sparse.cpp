#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
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
//' @title rcpp_jaccard_sparse
//' @name rcpp_jaccard_sparse
//' @description returns jaccard
//' @return list
//' @param kmer_counts_q list [mandatory]
//' @param kmer_counts_t list [mandatory]
//' @param k kmer length [default: 6]
//' @param min_jaccard min jaccard distance to report [default: 0.01]
//' @param sparse_threshold use kmer subset to evaluate change search strategy [default: 0.1]
//' @param sparse_n number of kmer subset to get sparse threshold [default: 20]
//' @param ncores number of cores [default: 1]
//' @param debug print to console [default: FALSE]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' l <- hiv |>
//'     MSA2dist::cds2aa() |>
//'     korthoR::count_kmers(k=6)
//' d <- korthoR::rcpp_jaccard_sparse(
//'     kmer_counts_q=l,
//'     kmer_counts_t=l,
//'     k=6,
//'     min_jaccard=0.01,
//'     sparse_threshold=0.1,
//'     sparse_n=100,
//'     ncores=1,
//'     debug=FALSE)
//' d
//' @export rcpp_jaccard_sparse
//' @author Kristian K Ullrich

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_jaccard_sparse(
  Rcpp::List kmer_counts_q,
  Rcpp::List kmer_counts_t,
  int k=6,
  double min_jaccard=0.01,
  double sparse_threshold=0.1,
  int sparse_n=20,
  int ncores=1,
  bool debug=false) {
  //get sizes
  int n = kmer_counts_q.size();
  int m = kmer_counts_t.size();
  if (debug) {
    std::cout << "number of proteins Q: " << n;
    std::cout << std::endl;
    std::cout << "number of proteins T: " << m;
    std::cout << std::endl;
  }
  //create std::vector for names, kmer counts, kmers
  //names
  std::vector<std::string> seqnames_q = kmer_counts_q.attr("names");
  std::vector<std::string> seqnames_t = kmer_counts_t.attr("names");
  //kmers
  std::vector<std::vector<std::string>> Qkmers_sorted(n);
  std::vector<std::vector<std::string>> Tkmers_sorted(m);
  //kmer counts
  std::vector<std::vector<int>> Qkmers_counts_sorted(n);
  std::vector<std::vector<int>> Tkmers_counts_sorted(m);
  //convert Rcpp vectors into std::vector
  Qkmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  Tkmers_sorted = convertInnerNamesRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  Qkmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_q);
  Tkmers_counts_sorted = convertRcppListOfIntegerVectorToStdVector(kmer_counts_t);
  //create kmerMaps
  std::map<std::string, std::vector<int>> QkmerMap;
  std::map<std::string, std::vector<int>> TkmerMap;
  //create vector to keep separate maps
  auto start_QkmerMap = std::chrono::steady_clock::now();
  std::vector<std::map<std::string, int>> QkmerMap_n(n);
  QkmerMap_n.resize(n);
  RcppThread::parallelFor(0, n, [&] (int i) {
    std::map<std::string, int> QkmerMap_i;
    std::vector<std::string> Qkmers_i_sorted = Qkmers_sorted[i];
    QkmerMap_i = vectorToMap(Qkmers_i_sorted, i);
    QkmerMap_n[i] = QkmerMap_i;
  }, ncores);
  //combine maps
  QkmerMap = combineMaps(QkmerMap_n);
  auto end_QkmerMap = std::chrono::steady_clock::now();
  auto duration_QkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_QkmerMap - start_QkmerMap);
  if (debug) {
    std::cout << "Time taken: QkmerMap creation " << duration_QkmerMap.count() << " milliseconds" << std::endl;
    std::cout << "QkmerMap_n size: " << QkmerMap_n.size() << std::endl;
    std::cout << "number of kmers QkmerMap: " << QkmerMap.size() << std::endl;
  }
  //create vector to keep separate maps
  auto start_TkmerMap = std::chrono::steady_clock::now();
  std::vector<std::map<std::string, int>> TkmerMap_m(m);
  TkmerMap_m.resize(m);
  RcppThread::parallelFor(0, m, [&] (int j) {
    std::map<std::string, int> TkmerMap_j;
    std::vector<std::string> Tkmers_j_sorted = Tkmers_sorted[j];
    TkmerMap_j = vectorToMap(Tkmers_j_sorted, j);
    TkmerMap_m[j] = TkmerMap_j;
  }, ncores);
  //combine maps
  TkmerMap = combineMaps(TkmerMap_m);
  auto end_TkmerMap = std::chrono::steady_clock::now();
  auto duration_TkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_TkmerMap - start_TkmerMap);
  if (debug) {
    std::cout << "Time taken: TkmerMap creation " << duration_TkmerMap.count() << " milliseconds" << std::endl;
    std::cout << "TkmerMap_n size: " << TkmerMap_m.size() << std::endl;
    std::cout << "number of kmers TkmerMap: " << TkmerMap.size() << std::endl;
  }
  //check sparse
  bool use_sparse = false;
  auto start_sparseCheck = std::chrono::steady_clock::now();
  std::map<std::string, std::vector<int>> subsetQkmerMap;
  subsetQkmerMap = getSubsetMap(QkmerMap, sparse_n);
  std::map<std::string, std::vector<int>> subsetTkmerMap;
  subsetTkmerMap = getSubsetMap(TkmerMap, sparse_n);
  std::vector<std::vector<int>> sparse_pairs(n);
  sparse_pairs = findPairs_sparse(subsetQkmerMap, subsetTkmerMap, n);
  std::tuple<std::vector<int>, std::vector<int>> sparsePairs;
  sparsePairs = flatten2tuple(sparse_pairs);
  std::vector<int> sparsePairs_i = std::get<0>(sparsePairs);
  std::vector<int> sparsePairs_j = std::get<1>(sparsePairs);
  double sparseValue = (static_cast<double>(sparsePairs_i.size())/(n*m));
  if (sparseValue<sparse_threshold) {
    use_sparse = true;
  }
  if (debug) {
    std::cout << "number of sparse candidate pairs: " << sparsePairs_i.size();
    std::cout << std::endl;
    if (use_sparse) {
      std::cout << "sparse value: " << sparseValue << " < sparse threshold: " << sparse_threshold << " >>> search strategy one vs one";
      std::cout << std::endl;
    } else {
      std::cout << "sparse value: " << sparseValue << " > sparse threshold: " << sparse_threshold << " >>> search strategy one vs many";
      std::cout << std::endl;
    }
  }
  auto end_sparseCheck = std::chrono::steady_clock::now();
  auto duration_sparseCheck = std::chrono::duration_cast<std::chrono::milliseconds>(end_sparseCheck - start_sparseCheck);
  if (debug) {
    std::cout << "Time taken: check sparse threshold " << duration_sparseCheck.count() << " milliseconds" << std::endl;
  }
  //get candidate pairs
  auto start_getCandidates = std::chrono::steady_clock::now();
  std::vector<std::vector<int>> Q_T_pairs(n);
  if (QkmerMap.size() < TkmerMap.size()) {
    if (use_sparse) {
      Q_T_pairs = findPairs_sparse(QkmerMap, TkmerMap, n);
    } else {
      RcppThread::parallelFor(0, n, [&] (int i) {
        Q_T_pairs[i] = findPairs(QkmerMap_n[i], TkmerMap);
      }, ncores);
    }
  } else {
    std::vector<std::vector<int>> T_Q_pairs(m);
    if (use_sparse) {
      T_Q_pairs = findPairs_sparse(TkmerMap, QkmerMap, m);
    } else {
      RcppThread::parallelFor(0, m, [&] (int j) {
        T_Q_pairs[j] = findPairs(TkmerMap_m[j], QkmerMap);
      }, ncores);
    }
    //transpose T_Q_hits to get Q_T_hits
    Q_T_pairs = transposeTQ(T_Q_pairs, n);
  }
  std::vector<std::set<int>> candidatesQT;
  candidatesQT = getCandidates(Q_T_pairs);
  std::vector<int> candidatePairs_ni(candidatesQT[0].begin(), candidatesQT[0].end());
  std::vector<int> candidatePairs_mj(candidatesQT[1].begin(), candidatesQT[1].end());
  int candidatePairs_ni_mj_size;
  candidatePairs_ni_mj_size = *candidatesQT[2].begin();
  auto end_getCandidates = std::chrono::steady_clock::now();
  auto duration_getCandidates = std::chrono::duration_cast<std::chrono::milliseconds>(end_getCandidates - start_getCandidates);
  if (debug) {
    std::cout << "Time taken: candidate pairs creation " << duration_getCandidates.count() << " milliseconds" << std::endl;
    std::cout << "number of candidate pairs: " << candidatePairs_ni_mj_size;
    std::cout << std::endl;
  }
  //calculate distances for each candidate pair
  auto start_calcDist = std::chrono::steady_clock::now();
  std::vector<std::vector<std::vector<double>>> Qni_Tmj_distances(
    candidatePairs_ni.size(),
    std::vector<std::vector<double>>(
      candidatePairs_mj.size(),
      std::vector<double>(
        6,
        0.0)));
  RcppThread::ThreadPool pool(ncores);
  pool.parallelFor(0, candidatePairs_ni.size(), [&] (int ni) {
    int Qni = candidatePairs_ni[ni];
    std::vector<int>Qni_T_pairs(Q_T_pairs[Qni]);
      pool.parallelFor(0, Qni_T_pairs.size(), [&, ni, Qni, Qni_T_pairs] (int mj) {
      int Tmj = Qni_T_pairs[mj];
      Qni_Tmj_distances[ni][mj]=getJaccardByIntegerVector(
        Qkmers_sorted[Qni],
        Tkmers_sorted[Tmj],
        Qkmers_counts_sorted[Qni],
        Tkmers_counts_sorted[Tmj],
        Qni,
        Tmj,
        k);
    });
  });
  pool.join();
  auto end_calcDist = std::chrono::steady_clock::now();
  auto duration_calcDist = std::chrono::duration_cast<std::chrono::milliseconds>(end_calcDist - start_calcDist);
  if (debug) {
    std::cout << "Time taken: distance calculation " << duration_calcDist.count() << " milliseconds" << std::endl;
  }
  std::vector<std::string> out_qname;
  std::vector<std::string> out_tname;
  std::vector<double> out_jaccard;
  std::vector<double> out_mash;
  std::vector<double> out_ani;
  std::vector<double> out_sumdist;
  for (const auto& row : Qni_Tmj_distances) {
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
