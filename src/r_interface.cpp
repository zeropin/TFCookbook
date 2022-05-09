// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "MOODS/moods.h"
#include "MOODS/moods_misc.h"
#include "MOODS/moods_tools.h"
#include "MOODS/match_types.h"
#include "MOODS/motif.h"
#include "MOODS/scanner.h"

using namespace Rcpp;

std::vector<std::vector<double> > mat_conversion(List& mats, size_t i){
  arma::mat tmp = as<arma::mat>(mats[i]);
  score_matrix out(tmp.n_rows);
    for (size_t i = 0; i < tmp.n_rows; ++i) {
        out[i] = arma::conv_to< std::vector<double> >::from(tmp.row(i));
    };
  return out;
}


std::vector<double>  get_thresholds(List mats,
          const std::vector<double> nuc_freqs,
          const double p){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[n + i] = thresholds[i];
  }
  return thresholds;
}


// [[Rcpp::export(.get_motif_positions)]]
List get_motif_positions(List& mats, const std::vector<std::string> x,
          const std::vector<double> nuc_freqs,
          const double E,
          const size_t w){
  size_t n = mats.size();
  std::vector<double> thresholds(2 * n);
  std::vector<score_matrix> matrices(2 * n);
  for(size_t i = 0; i < n; i++) {
    matrices[i] = mat_conversion(mats, i);
    matrices[n+i] = MOODS::tools::reverse_complement(matrices[i]);
    //thresholds[i] = MOODS::tools::threshold_from_p(matrices[i], nuc_freqs, p);
    thresholds[i] = E;
    thresholds[n + i] = thresholds[i];
  }
  MOODS::scan::Scanner scanner = MOODS::scan::Scanner(w);
  scanner.set_motifs(matrices, nuc_freqs, thresholds);
  size_t nstrings = x.size();
  std::vector<unsigned int> iloc;
  std::vector<unsigned int> jloc;
  std::vector<double> values;
  std::vector<char> strand;
  std::vector<size_t> pos;
  for (size_t i = 0; i < nstrings; i++){
    auto results = scanner.scan(x[i]);
    for (size_t j = 0; j < n; j++){
      if (results[j].size() > 0){
        for (size_t k = 0; k < results[j].size(); k++){
          iloc.push_back(i);
          jloc.push_back(j);
          values.push_back(results[j][k].score);
          strand.push_back('+');
          pos.push_back(results[j][k].pos);
        }
      }
      if (results[n+j].size() > 0){
        for (size_t k = 0; k < results[n+j].size(); k++){
          iloc.push_back(i);
          jloc.push_back(j);
          values.push_back(results[n+j][k].score);
          strand.push_back('-');
          pos.push_back(results[n+j][k].pos);
        }
      }
    }
  }
  return List::create(Named("motif_ix") = jloc,
                      Named("seq_ix") = iloc,
                      Named("strand") = strand,
                      Named("pos") = pos,
                      Named("score") = values);
}

