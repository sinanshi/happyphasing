#' @title An accurate and efficient statistical method for phasing haplotypes using large phased reference panels
#' @author Zhangyi He
#

#install.packages("plyr")
library("plyr")

source("HE2017_rfun_N0500_L0500.R")

load("h_chips.RData")
ref_pnl <- t(h)

args = commandArgs(trailingOnly = TRUE)
nrun <- as.integer(args[1])

N_pnl <- 1000
L_pnl <- 500

N_smp <- 20
L_smp <- 500

obs_dip <- array(NA, dim = c(N_smp, 2, L_smp))
obs_smp <- matrix(NA, nrow = N_smp, ncol = L_smp)
for (i in 1:N_smp) {
  obs_dip[i, , ] <- ref_pnl[(nrow(ref_pnl) - 2 * (i - 1) - 1):(nrow(ref_pnl) - 2 * (i - 1)), 1:L_smp]
  obs_smp[i, ] <- colSums(obs_dip[i, , ])
}

N_eff <- 15000
theta <- cmpgenerateSwitchRate_HRC_Chr20(ref_pnl, N_pnl, N_eff)[1:(L_pnl - 1)]

lambda <- cmpgenerateLambda(N_pnl)

ref_pnl <- ref_pnl[1:N_pnl, 1:L_pnl]

dip_smp_size <- 5

bck_len <- 50

if (nrun == 1) {
  best_dip_red_na <- cmprunBWPS_Diplo_Red_Pnl_NA(obs_smp, ref_pnl, bck_len, theta, lambda, dip_smp_size)
  save(best_dip_na, file = "dip_N0500_L0500_NA.RData")
} 
if (nrun == 2) {
  best_dip_red_fa <- cmprunBWPS_Diplo_Red_Pnl_FA(obs_smp, ref_pnl, bck_len, theta, lambda, dip_smp_size)
  save(best_dip_fa, file = "dip_N0500_L0500_FA.RData")
} 
if (nrun == 3) {
  best_dip_red_sa <- cmprunBWPS_Diplo_Red_Pnl_SA(obs_smp, ref_pnl, bck_len, theta, lambda, dip_smp_size)
  save(best_dip_sa, file = "dip_N0500_L0500_SA.RData")
} 
if (nrun == 4) { 
  best_dip_red_ta <- cmprunBWPS_Diplo_Red_Pnl_TA(obs_smp, ref_pnl, bck_len, theta, lambda, dip_smp_size)
  save(best_dip_ta, file = "dip_N0500_L0500_TA.RData")
}
