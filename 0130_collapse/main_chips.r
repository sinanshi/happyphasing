#' @title An accurate and efficient statistical method for phasing haplotypes using large phased reference panels
#' @author Zhangyi He
#
# something new: test
#install.packages("plyr")
library("plyr")

source("HE2017_rfun.R")
load("h_chips.RData")
ref_pnl <- t(h)

args = commandArgs(trailingOnly = TRUE)
nrun <- as.integer(args[1])
print(nrun)

N_pnl <- 2000
L_pnl <- 750
ref_pnl_N02000_L0750 <- ref_pnl[1:N_pnl, 1:L_pnl]
N_smp <- 1000
obs_smp_N02000_L0750 <- ref_pnl[(nrow(ref_pnl) - 20 + 1):nrow(ref_pnl), 1:L_pnl]

N_eff <- 15000
theta <- generateSwitchRate_HRC_Chr20(ref_pnl, N_pnl, N_eff)[1:(L_pnl - 1)]

lambda <- generateLambda(N_pnl)

block_length <- 50

M_sample <- 1


print(theta)
print(lambda)
#' Perform the forward procedure using the haplotype reference panel with state-space reduction
if (nrun == 1) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        print(obs_individual)
        dip_N02000_L0750[k, , ] <- performBWPS(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx0_set1.RData")
} 
if (nrun == 2) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx1(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx1_set1.RData")
} 
if (nrun == 3) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx2(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx2_set1.RData")
} 
if (nrun == 4) { 
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx3(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx3_set1.RData")
}


N_pnl <- 2000
L_pnl <- 750
ref_pnl_N02000_L0750 <- ref_pnl[1:N_pnl, 1:L_pnl]
N_smp <- 1000
obs_smp_N02000_L0750 <- ref_pnl[(nrow(ref_pnl) - 40 + 1):(nrow(ref_pnl) - 20), 1:L_pnl]

N_eff <- 15000
theta <- generateSwitchRate_HRC_Chr20(ref_pnl, N_pnl, N_eff)[1:(L_pnl - 1)]

lambda <- generateLambda(N_pnl)

block_length <- 50

M_sample <- 1

#' Perform the forward procedure using the haplotype reference panel with state-space reduction
if (nrun == 5) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx0_set2.RData")
} 
if (nrun == 6) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx1(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx1_set2.RData")
} 
if (nrun == 7) {
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx2(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx2_set2.RData")
} 
if (nrun == 8) { 
    dip_N02000_L0750 <- array(NA, dim = c(10, 2, L_pnl))
    for (k in 1:10) {
        obs_individual <- colSums(obs_smp_N02000_L0750[(2 * k - 1):(2 * k), ])
        dip_N02000_L0750[k, , ] <- performBWPS_Approx3(obs_individual, ref_pnl_N02000_L0750, theta, lambda, block_length, M_sample)
    }
    save(dip_N02000_L0750, file = "dip_N02000_L0750_Approx3_set2.RData")
}

