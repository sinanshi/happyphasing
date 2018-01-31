#' @title An accurate and efficient statistical method for phasing haplotypes using large phased reference panels
#' @author Zhangyi He
#

#install.packages("plyr")
library("plyr")

source("HE2017_rfun.R")
load("h_64940.RData")
ref_pnl <- t(h)

args = commandArgs(trailingOnly=TRUE)
nrun <- as.integer(args[1])
print(nrun)


#' A reference pnael created with 1000 individuals of 1000 markers (2.32 Mb) from the HRC (Chromosome 20)
N_pnl <- 2000
L_pnl <- 1000
ref_pnl_N02000_L1000 <- ref_pnl[1:N_pnl, 1:L_pnl]
N_smp <- 1000
obs_smp_N02000_L1000 <- ref_pnl[(nrow(ref_pnl) - N_smp * 2 + 1):nrow(ref_pnl), 1:L_pnl]

N_eff <- 15000
theta <- generateSwitchRate_HRC_Chr20(ref_pnl, N_pnl, N_eff)[1:(L_pnl - 1)]

lambda <- generateLambda(N_pnl)

block_length <- 50

M_sample <- 1

#' Perform the forward procedure using the haplotype reference panel with state-space reduction
obs_individual <- colSums(obs_smp_N02000_L1000[1:2, ])

if (nrun == 1) { 
    dip_N02000_L1000 <- performBWPS(obs_individual, ref_pnl_N02000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N02000_L1000, file="dip_N02000_L1000.RData")
} 

if (nrun == 2) { 
    dip_N02000_L1000_Approx1 <- performBWPS_Approx1(obs_individual, ref_pnl_N02000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N02000_L1000_Approx1, file="dip_N02000_L1000_Approx1.RData")
} 
if (nrun == 3) { 
    dip_N02000_L1000_Approx2 <- performBWPS_Approx2(obs_individual, ref_pnl_N02000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N02000_L1000_Approx2, file="dip_N02000_L1000_Approx2.RData")
} 

if (nrun == 4) { 
    dip_N02000_L1000_Approx3 <- performBWPS_Approx3(obs_individual, ref_pnl_N02000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N02000_L1000_Approx3, file="dip_N02000_L1000_Approx3.RData")
}
#' A reference pnael created with 2000 individuals of 1000 markers (2.32 Mb) from the HRC (Chromosome 20)
N_pnl <- 4000
L_pnl <- 1000
ref_pnl_N04000_L1000 <- ref_pnl[1:N_pnl, 1:L_pnl]
N_smp <- 1000
obs_smp_N04000_L1000 <- ref_pnl[(nrow(ref_pnl) - N_smp * 2 + 1):nrow(ref_pnl), 1:L_pnl]

N_eff <- 15000
theta <- generateSwitchRate_HRC_Chr20(ref_pnl, N_pnl, N_eff)[1:(L_pnl - 1)]

lambda <- generateLambda(N_pnl)

block_length <- 50

M_sample <- 1

#' Perform the forward procedure using the haplotype reference panel with state-space reduction
obs_individual <- colSums(obs_smp_N04000_L1000[1:2, ])

if (nrun == 5) { 
    dip_N04000_L1000 <- performBWPS(obs_individual, ref_pnl_N04000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N04000_L1000, file="dip_N04000_L1000.RData")
} 

if (nrun == 6) { 
    dip_N04000_L1000_Approx1 <- performBWPS_Approx1(obs_individual, ref_pnl_N04000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N04000_L1000_Approx1, file="dip_N04000_L1000_Approx1.RData")
} 

if (nrun == 7) { 
    dip_N04000_L1000_Approx2 <- performBWPS_Approx2(obs_individual, ref_pnl_N04000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N04000_L1000_Approx2, file="dip_N04000_L1000_Approx2.RData")
} 

if (nrun == 8) { dip_N04000_L1000_Approx3 <- performBWPS_Approx3(obs_individual, ref_pnl_N04000_L1000, theta, lambda, block_length, M_sample)
    save(dip_N04000_L1000_Approx3, file="dip_N04000_L1000_Approx3.RData")
} 

