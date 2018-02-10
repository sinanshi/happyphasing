// A fast and accurate statistical method for phasing haplotypes using large phased reference panels
// Zhangyi He
//

// C functions

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;
using namespace std;

// Calculate the emission probabilities at a specific marker of the genome
// [[Rcpp::export]]
NumericMatrix calculateEmnProb_Diplo_cpp(int obs_gen, IntegerVector ref_hap_a, IntegerVector ref_hap_b, NumericMatrix emn_prob_tab) {
  int N_ref_hap_a = ref_hap_a.length();
  int N_ref_hap_b = ref_hap_b.length();
  
  int ref_gen;
  
  NumericMatrix emn_prob(N_ref_hap_a, N_ref_hap_b);
  for (int i = 0; i < N_ref_hap_a; i++) {
    for (int j = 0; j < N_ref_hap_b; j++) {
      ref_gen = ref_hap_a(i) + ref_hap_b(j);
      emn_prob(i, j) = emn_prob_tab(ref_gen, obs_gen);
    }
  }
  
  return(emn_prob);
}


// Collapse the collapsed forward probabilities at the right boundary of the genomic block for case 1
// [[Rcpp::export]]
NumericMatrix collapseFwdProb_Diplo_FullRed_Bck_NA_cpp(NumericMatrix uncollapsed_fwd_prob_lbndry, int N_pnl, int N_bck, IntegerVector collapsed_map) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix collapsed_fwd_prob_lbndry(N_bck, N_bck);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int v = collapsed_map(j) - 1;
      collapsed_fwd_prob_lbndry(u, v) = collapsed_fwd_prob_lbndry(u, v) + uncollapsed_fwd_prob_lbndry(i, j);
    }
  }
  
  return(collapsed_fwd_prob_lbndry);
}

// Uncollapse the collapsed forward probabilities at the right boundary of the genomic block for case 1
// [[Rcpp::export]]
NumericMatrix uncollapseFwdProb_Diplo_FullRed_Bck_NA_cpp(NumericMatrix collapsed_fwd_prob_nr_nr_rbndry, NumericMatrix collapsed_fwd_prob_re_re_rbndry, 
                                                         NumericMatrix collapsed_fwd_prob_lbndry, NumericMatrix uncollapsed_fwd_prob_lbndry, 
                                                         int N_pnl, NumericVector collapsed_map, NumericVector collapsed_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix uncollapsed_fwd_prob_rbndry(N_pnl, N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int cnt_bck_hap_a = collapsed_cnt(u);
      int v = collapsed_map(j) - 1;
      int cnt_bck_hap_b = collapsed_cnt(v);
      uncollapsed_fwd_prob_rbndry(i, j) = 
        collapsed_fwd_prob_nr_nr_rbndry(u, v) * uncollapsed_fwd_prob_lbndry(i, j) / collapsed_fwd_prob_lbndry(u, v) + 
        collapsed_fwd_prob_re_re_rbndry(u, v) / cnt_bck_hap_a / cnt_bck_hap_b;
    }
  }
  
  return(uncollapsed_fwd_prob_rbndry);
}

// Collapse the collapsed forward probabilities at the right boundary of the genomic block for case 1
// [[Rcpp::export]]
NumericMatrix collapseFwdProb_Diplo_SemiRed_Bck_NA_cpp(NumericMatrix uncollapsed_fwd_prob_lbndry, int N_pnl, int N_bck, IntegerVector collapsed_map) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix collapsed_fwd_prob_lbndry(N_pnl, N_bck);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int v = collapsed_map(j) - 1;
      collapsed_fwd_prob_lbndry(i, v) = collapsed_fwd_prob_lbndry(i, v) + uncollapsed_fwd_prob_lbndry(i, j);
    }
  }
  
  return(collapsed_fwd_prob_lbndry);
}

// Uncollapse the collapsed forward probabilities at the right boundary of the genomic block for case 2
// [[Rcpp::export]]
NumericMatrix uncollapseFwdProb_Diplo_SemiRed_Bck_NA_cpp(NumericMatrix collapsed_fwd_prob_nr_re_rbndry, 
                                                         int N_pnl, IntegerVector collapsed_map, IntegerVector collapsed_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix uncollapsed_fwd_prob_rbndry(N_pnl, N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int cnt_bck_hap_a = collapsed_cnt(u);
      int v = collapsed_map(j) - 1;
      int cnt_bck_hap_b = collapsed_cnt(v);
      uncollapsed_fwd_prob_rbndry(i, j) = collapsed_fwd_prob_nr_re_rbndry(i, v) / cnt_bck_hap_b + collapsed_fwd_prob_nr_re_rbndry(j, u) / cnt_bck_hap_a;
    }
  }
  
  return(uncollapsed_fwd_prob_rbndry);
}


// collapse the uncollapsed forward probabilities at the right boundary of the genomic block using the 1st type of the approximation
// [[Rcpp::export]]
NumericMatrix collapseFwdProb_Diplo_FullRed_Bck_FA_cpp(NumericMatrix uncollapsed_fwd_prob_lbndry, int N_pnl, int N_bck, IntegerVector collapsed_map) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix collapsed_fwd_prob_lbndry(N_bck, N_bck);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int v = collapsed_map(j) - 1;
      collapsed_fwd_prob_lbndry(u, v) = collapsed_fwd_prob_lbndry(u, v) + uncollapsed_fwd_prob_lbndry(i, j);
    }
  }
  
  return(collapsed_fwd_prob_lbndry);
}

// Uncollapse the collapsed forward probabilities at the right boundary of the genomic block using the 1st type of the approximation
// [[Rcpp::export]]
NumericMatrix uncollapseFwdProb_Diplo_FullRed_Bck_FA_cpp(NumericMatrix collapsed_fwd_prob_nr_nr_rbndry, NumericMatrix collapsed_fwd_prob_re_re_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_nr_re_rbndry, NumericMatrix collapsed_fwd_prob_re_nr_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_lbndry, NumericMatrix uncollapsed_fwd_prob_lbndry,
                                                         int N_pnl, IntegerVector collapsed_map, IntegerVector collapsed_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  int N_bck = collapsed_fwd_prob_lbndry.nrow();
  
  NumericVector marginal_collapsed_fwd_prob_lbndry(N_bck);
  for (int u = 0; u < N_bck; u++) {
    marginal_collapsed_fwd_prob_lbndry(u) = sum(collapsed_fwd_prob_lbndry(u, _));
  }
  
  NumericVector marginal_uncollapsed_fwd_prob_lbndry(N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    marginal_uncollapsed_fwd_prob_lbndry(i) = sum(uncollapsed_fwd_prob_lbndry(i, _));
  }
  
  NumericMatrix uncollapsed_fwd_prob_rbndry(N_pnl, N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int cnt_bck_hap_a = collapsed_cnt(u);
      int v = collapsed_map(j) - 1;
      int cnt_bck_hap_b = collapsed_cnt(v);
      uncollapsed_fwd_prob_rbndry(i, j) =
        collapsed_fwd_prob_nr_nr_rbndry(u, v) * uncollapsed_fwd_prob_lbndry(i, j) / collapsed_fwd_prob_lbndry(u, v) +
        collapsed_fwd_prob_nr_re_rbndry(u, v) * (marginal_uncollapsed_fwd_prob_lbndry(i) / marginal_collapsed_fwd_prob_lbndry(u)) / cnt_bck_hap_b +
        collapsed_fwd_prob_re_nr_rbndry(u, v) / cnt_bck_hap_a * (marginal_uncollapsed_fwd_prob_lbndry(j) / marginal_collapsed_fwd_prob_lbndry(v)) +
        collapsed_fwd_prob_re_re_rbndry(u, v) / cnt_bck_hap_a / cnt_bck_hap_b;
    }
  }

  return(uncollapsed_fwd_prob_rbndry);
}


// collapse the uncollapsed forward probabilities at the right boundary of the genomic block using the 2nd type of the approximation
// [[Rcpp::export]]
NumericMatrix collapseFwdProb_Diplo_FullRed_Bck_SA_cpp(NumericMatrix uncollapsed_fwd_prob_lbndry, int N_pnl, int N_bck, IntegerVector collapsed_map) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix collapsed_fwd_prob_lbndry(N_bck, N_bck);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int v = collapsed_map(j) - 1;
      collapsed_fwd_prob_lbndry(u, v) = collapsed_fwd_prob_lbndry(u, v) + uncollapsed_fwd_prob_lbndry(i, j);
    }
  }
  
  return(collapsed_fwd_prob_lbndry);
}

// Uncollapse the collapsed forward probabilities at the right boundary of the genomic block using the 2nd type of the approximation
// [[Rcpp::export]]
NumericMatrix uncollapseFwdProb_Diplo_FullRed_Bck_SA_cpp(NumericMatrix collapsed_fwd_prob_nr_nr_rbndry, NumericMatrix collapsed_fwd_prob_re_re_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_nr_re_rbndry, NumericMatrix collapsed_fwd_prob_re_nr_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_lbndry, NumericMatrix uncollapsed_fwd_prob_lbndry,
                                                         int N_pnl, IntegerVector collapsed_map, IntegerVector collapsed_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  int N_bck = collapsed_fwd_prob_lbndry.nrow();
  
  NumericVector marginal_collapsed_fwd_prob_lbndry(N_bck);
  for (int u = 0; u < N_bck; u++) {
    marginal_collapsed_fwd_prob_lbndry(u) = sum(collapsed_fwd_prob_lbndry(u, _));
  }
  
  NumericVector marginal_uncollapsed_fwd_prob_lbndry(N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    marginal_uncollapsed_fwd_prob_lbndry(i) = sum(uncollapsed_fwd_prob_lbndry(i, _));
  }
  
  NumericMatrix uncollapsed_fwd_prob_rbndry(N_pnl, N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    for (int j = 0; j < N_pnl; j++) {
      int u = collapsed_map(i) - 1;
      int cnt_bck_hap_a = collapsed_cnt(u);
      int v = collapsed_map(j) - 1;
      int cnt_bck_hap_b = collapsed_cnt(v);
      uncollapsed_fwd_prob_rbndry(i, j) =
        collapsed_fwd_prob_nr_nr_rbndry(u, v) * (marginal_uncollapsed_fwd_prob_lbndry(i) / marginal_collapsed_fwd_prob_lbndry(u)) * (marginal_uncollapsed_fwd_prob_lbndry(j) / marginal_collapsed_fwd_prob_lbndry(v)) +
        collapsed_fwd_prob_nr_re_rbndry(u, v) * (marginal_uncollapsed_fwd_prob_lbndry(i) / marginal_collapsed_fwd_prob_lbndry(u)) / cnt_bck_hap_b +
        collapsed_fwd_prob_re_nr_rbndry(u, v) / cnt_bck_hap_a * (marginal_uncollapsed_fwd_prob_lbndry(j) / marginal_collapsed_fwd_prob_lbndry(v)) +
        collapsed_fwd_prob_re_re_rbndry(u, v) / cnt_bck_hap_a / cnt_bck_hap_b;
    }
  }
  
  return(uncollapsed_fwd_prob_rbndry);
}


// collapse the uncollapsed forward probabilities at the right boundary of the genomic block using the 3rd type of the approximation
// [[Rcpp::export]]
NumericMatrix collapseFwdProb_Diplo_FullRed_Bck_TA_cpp(NumericVector marginal_uncollapsed_fwd_prob_lbndry, int N_pnl, int N_bck, IntegerVector collapsed_map) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericVector marginal_collapsed_fwd_prob_lbndry(N_bck);
  for (int i = 0; i < N_pnl; i++) {
    int u = collapsed_map(i) - 1;
    marginal_collapsed_fwd_prob_lbndry(u) = marginal_collapsed_fwd_prob_lbndry(u) + marginal_uncollapsed_fwd_prob_lbndry(i);
  }
  
  NumericMatrix collapsed_fwd_prob_lbndry(N_bck, N_bck);
  for (int u = 0; u < N_bck; u++) {
    for (int v = 0; v < N_bck; v++) {
      collapsed_fwd_prob_lbndry(u, v) = marginal_collapsed_fwd_prob_lbndry(u) * marginal_collapsed_fwd_prob_lbndry(v);
    }
  }
  
  return(collapsed_fwd_prob_lbndry);
}

// Uncollapse the collapsed forward probabilities at the right boundary of the genomic block using the 3rd type of the approximation
// [[Rcpp::export]]
NumericVector uncollapseFwdProb_Diplo_FullRed_Bck_TA_cpp(NumericMatrix collapsed_fwd_prob_nr_nr_rbndry, NumericMatrix collapsed_fwd_prob_re_re_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_nr_re_rbndry, NumericMatrix collapsed_fwd_prob_re_nr_rbndry,
                                                         NumericMatrix collapsed_fwd_prob_lbndry, NumericVector marginal_uncollapsed_fwd_prob_lbndry,
                                                         int N_pnl, IntegerVector collapsed_map, IntegerVector collapsed_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  int N_bck = collapsed_fwd_prob_lbndry.nrow();
  
  NumericVector marginal_collapsed_fwd_prob_lbndry(N_bck);
  
  NumericVector marginal_collapsed_fwd_prob_rbndry_nr(N_bck);
  NumericVector marginal_collapsed_fwd_prob_rbndry_re(N_bck);
  
  for (int u = 0; u < N_bck; u++) {
    marginal_collapsed_fwd_prob_lbndry(u) = sum(collapsed_fwd_prob_lbndry(u, _));
    
    marginal_collapsed_fwd_prob_rbndry_nr(u) = sum(collapsed_fwd_prob_nr_nr_rbndry(u, _) + collapsed_fwd_prob_nr_re_rbndry(u, _));
    marginal_collapsed_fwd_prob_rbndry_re(u) = sum(collapsed_fwd_prob_re_nr_rbndry(u, _) + collapsed_fwd_prob_re_re_rbndry(u, _));
  }
  
  NumericVector uncollapsed_fwd_prob_rbndry(N_pnl);
  for (int i = 0; i < N_pnl; i++) {
    int u = collapsed_map(i) - 1;
    int cnt_bck_hap_a = collapsed_cnt(u);
    uncollapsed_fwd_prob_rbndry(i) =
      marginal_collapsed_fwd_prob_rbndry_nr(u) * (marginal_uncollapsed_fwd_prob_lbndry(i) / marginal_collapsed_fwd_prob_lbndry(u)) +
      marginal_collapsed_fwd_prob_rbndry_re(u) / cnt_bck_hap_a;
  }
  
  return(uncollapsed_fwd_prob_rbndry);
}

