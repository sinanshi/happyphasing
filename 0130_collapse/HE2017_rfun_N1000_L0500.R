#' @title An accurate and efficient statistical method for phasing haplotypes using large phased reference panels
#' @author Zhangyi He
#'

#' R functions

#install.packages("data.table")
library("data.table")
#install.packages("abind")
library("abind")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppEigen")
library("RcppEigen")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

#' call C++ functions
sourceCpp("HE2017_cfun.cpp")

#############################################################################################################################

#' @section Section 1. Introduction

#############################################################################################################################

#' @section Section 2. Materials and Methods

#' Initialisation

#' Generate a haplotype reference panel
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param L_pnl the number of the markers across the chromosome in the reference panel
#' @return a haplotype reference panel returned in a matrix (row = haplotype and col = marker)

#' Standard version
generateRefPanel <- function(N_pnl, L_pnl, seed = 1) {
  set.seed(seed)
  
  ref_pnl <- matrix(as.integer(sample(c(0, 1), N_pnl * L_pnl, TRUE, rep(0.5, 2))), nrow = N_pnl, ncol = L_pnl)  
  
  return(ref_pnl)
}
#' Compiled version
cmpgenerateRefPanel <- cmpfun(generateRefPanel)

#########################

#' Generate the switch rates across the chromosome
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param N_eff the effective population size
#' @param d the genetic distances between successive markers
#' @return the switch rates across the chromosome returned in a vector

#' Standard version
generateSwitchRate <- function(N_pnl, N_eff, d) {
  rho <- 4 * N_eff * d
  theta <- 1 - exp(-rho / N_pnl)
  
  return(theta)
}
#' Compiled version
cmpgenerateSwitchRate <- cmpfun(generateSwitchRate)

#' Generate the switch rates across the chromosome 20 with the HRC
#' @param ref_pnl the reference panel for chromosome 20 created by the HRC
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param N_eff the effective population size
#' @return the switch rates across the chromosome 20 with the HRS returned in a vector

#' Standard version
generateSwitchRate_HRC_Chr20 <- function(ref_pnl, N_pnl, N_eff) {
  genetic_map <- read.table("genetic_map_chr20_combined_b37.txt", header = TRUE)
  genetic_map <- approx(as.vector(genetic_map[, 1]), as.vector(genetic_map[, 3]), xout = as.numeric(as.vector(colnames(ref_pnl))))
  genetic_map_cm <- genetic_map$y
  genetic_map_cm[which(genetic_map_cm < 0)] <- 0
  d <- diff(genetic_map_cm, lag = 1) / 100
  rho <- 4 * N_eff * d
  theta <- 1 - exp(-rho / N_pnl)
  
  return(theta)
}
#' Compiled version
cmpgenerateSwitchRate_HRC_Chr20 <- cmpfun(generateSwitchRate_HRC_Chr20)

#########################

#' Generate the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @return the mutation rate

#' Standard version
generateLambda <- function(N_pnl) {
  theta <- 1 / sum(1 / (1:(N_pnl - 1)))
  lambda <- 0.5 * theta / (N_pnl + theta)
  
  return(lambda)
}
#' Compiled version
cmpgenerateLambda <- cmpfun(generateLambda)

#########################

#' Generate a observed genotype sample
#' @param N_smp the number of the individuals in the observed sample
#' @param L_smp the number of the markers across the chromosome in the observed sample (= L_pnl)
#' @return a observed genotype sample returned in a matrix (row = genotype and col = marker)
#' N.B. The observed genotype could be one of the values from the set {0, 1, 2, 3}, where 0 = rr, 1 = ra/ar, 2 = aa, 3 = ?

#' Standard version
generateObsSample <- function(N_smp, L_smp, seed = 2) {
  set.seed(seed)
  
  obs_smp <- matrix(as.integer(sample(c(0, 1, 2, 3), N_smp * L_smp, TRUE, rep(0.25, 4))), nrow = N_smp, ncol = L_smp)
  
  return(obs_smp)
}
#' Compiled version
cmpgenerateObsSample <- cmpfun(generateObsSample)

###########################################################################

#' Haplotype phasing

#' Emission probability calculation

#' Generate the emission probability table using Tab.1 
#' @param lambda the mutation rate
#' @return the emission probability table returned in a matrix (row = hidden genotype, col = observed genotype)

#' Standard version
generateEmnProbTab <- function(lambda) {
  emn_prob <- matrix(NA, nrow = 3, ncol = 4)
  
  emn_prob[1, 1] <- (1 - lambda) ^ 2
  emn_prob[2, 1] <- lambda * (1 - lambda)
  emn_prob[3, 1] <- lambda ^ 2
  
  emn_prob[1, 2] <- 2 * lambda * (1 - lambda)
  emn_prob[2, 2] <- lambda ^ 2 + (1 - lambda) ^ 2
  emn_prob[3, 2] <- 2 * lambda * (1 - lambda)
  
  emn_prob[1, 3] <- lambda ^ 2
  emn_prob[2, 3] <- lambda * (1 - lambda)
  emn_prob[3, 3] <- (1 - lambda) ^ 2
  
  emn_prob[1, 4] <- 1
  emn_prob[2, 4] <- 1
  emn_prob[3, 4] <- 1
  
  return(emn_prob)
}
#' Compiled version
cmpgenerateEmnProbTab <- cmpfun(generateEmnProbTab)

#' Calculate the emission probabilities at a specific marker of the chromosome using Tab.1
#' @param obs_gen the observed genotype at a specific marker
#' @param ref_hap_a the 1st reference haplotypes at a specific marker
#' @param ref_hap_b the 2nd reference haplotypes at a specific marker
#' @param emn_prob_tab the emission probability table
#' @return the emission probabilities at a specific marker of the chromosome returned in a matrix (row = haplotype a, col = haplotype b)
#' N.B. The observed genotype could be one of the values from the set {0, 1, 2, 3}, where 0 = rr, 1 = ra/ar, 2 = aa, 3 = ?

#' Standard version
calculateEmnProb_Diplo <- function(obs_gen, ref_hap_a, ref_hap_b, emn_prob_tab) {
  emn_prob <- function(ref_gen) {
    return(emn_prob_tab[(ref_gen + 1), (obs_gen + 1)])
  }
  
  N_hap_a <- length(ref_hap_a)
  N_hap_b <- length(ref_hap_b)
  ref_gen <- matrix(ref_hap_a, nrow = N_hap_a, ncol = N_hap_b) + t(matrix(ref_hap_b, nrow = N_hap_b, ncol = N_hap_a))
  
  return(apply(ref_gen, 2, emn_prob))
}
#' Compiled version
cmpcalculateEmnProb_Diplo <- cmpfun(calculateEmnProb_Diplo)

#########################

#' Generate the sampled diplotype at a specific marker of the chromosome using Eqs.69-73
#' @param obs_gen the observed genotype at a specific marker
#' @param ref_dip the reference diplotype at a specific marker
#' @param lambda the mutation rate
#' @return the sampled diplotype at a specific marker of the chromosome returned in a vector

#' Standard version
generateDiplotype_Diplo <- function(obs_gen, ref_dip, lambda) {
  smp_dip <- numeric(2)
  
  if (obs_gen == 0) {
    smp_dip <- c(0, 0)
  }
  if (obs_gen == 1) {
    if (sum(ref_dip) == 1) {
      smp_prob <- c((1 - lambda) ^ 2, lambda ^ 2) / ((1 - lambda) ^ 2 + lambda ^ 2)
      choice <- sample.int(2, size = 1, replace = TRUE, prob = smp_prob)
      if (choice == 1) {
        smp_dip <- c(ref_dip[1], ref_dip[2])
      }
      if (choice == 2) {
        smp_dip <- c(1 - ref_dip[1], 1 - ref_dip[2])
      }
    }
    else {
      smp_prob <- c(0.5, 0.5)
      choice <- sample.int(2, size = 1, replace = TRUE, prob = smp_prob)
      if (choice == 1) {
        smp_dip <- c(0, 1)
      }
      if (choice == 2) {
        smp_dip <- c(1, 0)
      }
    }
  }
  if (obs_gen == 2) {
    smp_dip <- c(1, 1)
  }
  if (obs_gen == 3) {
    smp_prob <- c((1 - lambda) ^ 2, (1 - lambda) * lambda, lambda * (1 - lambda), lambda ^ 2)
    choice <- sample.int(4, size = 1, replace = TRUE, prob = smp_prob)
    if (choice == 1) {
      smp_dip <- c(ref_dip[1], ref_dip[2])
    }
    if (choice == 2) {
      smp_dip <- c(ref_dip[1], 1 - ref_dip[2])
    }
    if (choice == 3) {
      smp_dip <- c(1 - ref_dip[1], ref_dip[2])
    } 
    if (choice == 4) {
      smp_dip <- c(1 - ref_dip[1], 1 - ref_dip[2])
    }  
  }
  
  return(smp_dip)
}
#' Compiled version
cmpgenerateDiplotype_Diplo <- cmpfun(generateDiplotype_Diplo)

#########################

#' Construct the histogram of the sampled diplotypes across the chromosome using Eq.74
#' @param smp_dip the sampled diplotypes
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the histogram of the sampled diplotypes across the chromosome including the sampled diplotypes and the corresponding counts returned in a list

#' Standard version
constructHist_Diplo <- function(smp_dip, smp_dip_size) {
  idx <- NULL
  
  for (m in 1:smp_dip_size) {
    if (all(smp_dip[m, ,] == smp_dip[1, 1:2, ]) || all(smp_dip[m, , ] == smp_dip[1, 2:1, ])) {
      idx <- c(idx, m)
    }
  }
  dip <- array(smp_dip[1, , ], dim = c(1, dim(smp_dip)[2], dim(smp_dip)[3]))
  cnt <- length(idx)
  
  idx_set <- 1:smp_dip_size
  idx_set <- idx_set[!(idx_set %in% idx)]
  
  while (length(idx_set) > 0) {
    idx <- NULL
    
    for (m in idx_set) {
      if (all(smp_dip[m, ,] == smp_dip[idx_set[1], 1:2, ]) || all(smp_dip[m, , ] == smp_dip[idx_set[1], 2:1, ])) {
        idx <- c(idx, m)
      }
    }
    dip <- abind(dip, smp_dip[idx_set[1], , ], along = 1)
    cnt <- c(cnt, length(idx))
    
    idx_set <- idx_set[!(idx_set %in% idx)]
  }
  
  return(list(diplotype = dip, count = cnt))
}
#' Compiled version
cmpconstructHist_Diplo <- cmpfun(constructHist_Diplo)

##################################################

#' Haplotype phasing without state-space reduction

#' Calculate the transition probabilities at a specific marker of the chromosome (for uncollapsed state space) using Eq.4
#' @param ref_hap_a_idx the index of the reference haplotype a
#' @param ref_hap_b_idx the index of the reference haplotype b
#' @param theta the switch rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @return the transition probabilities at a specific marker of the chromosome returned in a matrix (row = haplotype a, col = haplotype b)

#' Standard version
calculateTransProb_Diplo_Std_Pnl <- function(ref_hap_a_idx, ref_hap_b_idx, theta, N_pnl) {
  trans_prob <- matrix(rep((theta / N_pnl) * (theta / N_pnl), N_pnl * N_pnl), nrow = N_pnl, ncol = N_pnl)
  trans_prob[ref_hap_a_idx, ] <- rep((theta / N_pnl + (1 - theta)) * (theta / N_pnl), N_pnl)
  trans_prob[, ref_hap_b_idx] <- rep((theta / N_pnl) * (theta / N_pnl + (1 - theta)), N_pnl)
  trans_prob[ref_hap_a_idx, ref_hap_b_idx] <- (theta / N_pnl + (1 - theta)) * (theta / N_pnl + (1 - theta))
  
  return(trans_prob)
}
#' Compiled version
cmpcalculateTransProb_Diplo_Std_Pnl <- cmpfun(calculateTransProb_Diplo_Std_Pnl)

#########################

#' Forward procedure without state-space reduction

#' Run the forward procedure without state-space reduction using Eqs.8-10
#' @param obs_gen the observed individual
#' @param ref_pnl the reference panel
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @return the forward probabilities across the chromosome returned in an array (dim = c(haplotype a, haplotype b, marker))

#' Standard version
runFwdProcedure_Diplo_Std_Pnl <- function(obs_gen, ref_pnl, theta, lambda) {
  message("forward procedure")
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  
  emn_prob_tab <- cmpgenerateEmnProbTab(lambda)
  
  fwd_prob <- array(NA, dim = c(N_pnl, N_pnl, L_pnl))
  
  # initialise the forward probabilities at the first marker of the chromosome using Eq.10
  # print(paste("marker:", 1))
  # fwd_prob[, , 1] <- cmpcalculateEmnProb_Diplo(obs_gen[1], ref_pnl[, 1], ref_pnl[, 1], emn_prob_tab) / (N_pnl ^ 2)
  fwd_prob[, , 1] <- calculateEmnProb_Diplo_cpp(obs_gen[1], ref_pnl[, 1], ref_pnl[, 1], emn_prob_tab) / (N_pnl ^ 2)
  
  # update the forward probabilities across the chromosome (from the first marker to the last marker) using Eq.9
  if (L_pnl > 1) {
    for (l in 2:L_pnl) {
      # print(paste("marker:", l))
      # emn_prob <- cmpcalculateEmnProb_Diplo(obs_gen[l], ref_pnl[, l], ref_pnl[, l], emn_prob_tab)
      emn_prob <- calculateEmnProb_Diplo_cpp(obs_gen[l], ref_pnl[, l], ref_pnl[, l], emn_prob_tab)
      sgl_sum_hap <- rowSums(fwd_prob[, , l - 1])
      dbl_sum_hap <- sum(fwd_prob[, , l - 1])
      ns_rate <- 1 - theta[l - 1]
      sw_rate <- theta[l - 1] / N_pnl
      for (i in 1:N_pnl) {
        for (j in 1:N_pnl) {
          fwd_prob[i, j, l] <- emn_prob[i, j] * 
            ((ns_rate ^ 2) * fwd_prob[i, j, l - 1] + ns_rate * sw_rate * sgl_sum_hap[i] + sw_rate * ns_rate * sgl_sum_hap[j] + (sw_rate ^ 2) * dbl_sum_hap)
        }
      }
      if (min(fwd_prob[, , l]) < 1e-30) {
        fwd_prob[, , l] <- fwd_prob[, , l] * 1e+20
      }
      if (max(fwd_prob[, , l]) > 1e+30) {
        fwd_prob[, , l] <- fwd_prob[, , l] * 1e-20
      }
    }
  }
  
  return(fwd_prob)
}
#' Compiled version
cmprunFwdProcedure_Diplo_Std_Pnl <- cmpfun(runFwdProcedure_Diplo_Std_Pnl)

#########################

#' Backward procedure without state-space reduction

#' Run the backward procedure without state-space reduction using Eqs.60-61 and 68-73
#' @param obs_gen the observed individual
#' @param ref_pnl the reference panel
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param fwd_prob the forward probabilities across the chromosome
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the sampled diplotypes across the chromosome returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBkdProcedure_Diplo_Std_Pnl <- function(obs_gen, ref_pnl, theta, lambda, fwd_prob, smp_dip_size) {
  message("backward procedure")
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  
  bkd_path <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  
  # initialise the sampled diplotypes at the last marker of the chromosome using Eqs.60-61
  # print(paste("marker:", L_pnl))
  for (m in 1:smp_dip_size) {
    smp_prob <- as.vector(fwd_prob[, , L_pnl] / sum(fwd_prob[, , L_pnl]))
    bkd_path_idx <- sample.int(N_pnl * N_pnl, size = 1, replace = TRUE, prob = smp_prob)
    bkd_path[m, , L_pnl] <- if (bkd_path_idx %% N_pnl == 0) c(N_pnl, bkd_path_idx %/% N_pnl) else c(bkd_path_idx %% N_pnl, bkd_path_idx %/% N_pnl + 1)
    smp_dip[m, , L_pnl] <- cmpgenerateDiplotype_Diplo(obs_gen[L_pnl], ref_pnl[bkd_path[m, , L_pnl], L_pnl], lambda)
  }
  
  # update the sampled diplotypes across the chromosome (from the last marker to the first marker) using Eqs.60-61
  if (L_pnl > 1) {
    for (l in (L_pnl - 1):1) {
      # print(paste("marker:", l))
      for (m in 1:smp_dip_size) {
        trans_prob <- cmpcalculateTransProb_Diplo_Std_Pnl(bkd_path[m, 1, l + 1], bkd_path[m, 2, l + 1], theta[l], N_pnl)
        smp_prob <- as.vector(fwd_prob[, , l] * trans_prob / sum(fwd_prob[, , l] * trans_prob))
        bkd_path_idx <- sample.int(N_pnl * N_pnl, size = 1, replace = TRUE, prob = smp_prob)
        bkd_path[m, , l] <- if (bkd_path_idx %% N_pnl == 0) c(N_pnl, bkd_path_idx %/% N_pnl) else c(bkd_path_idx %% N_pnl, bkd_path_idx %/% N_pnl + 1)
        smp_dip[m, , l] <- cmpgenerateDiplotype_Diplo(obs_gen[l], ref_pnl[bkd_path[m, , l], l], lambda)
      }
    }
  }
  
  #return(list(path = bkd_path, diplotype = smp_dip))
  return(smp_dip)
}
#' Compiled version
cmprunBkdProcedure_Diplo_Std_Pnl <- cmpfun(runBkdProcedure_Diplo_Std_Pnl)

#########################

#' Run the block-wise path sampling without state-space reduction using Eqs.8-10, 60-61 and 68-74
#' @param obs_smp the observed sample
#' @param ref_pnl the reference panel
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the estimated diplotypes across the chromosome achieving the maximum a posterior returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBWPS_Diplo_Std_Pnl <- function(obs_smp, ref_pnl, theta, lambda, smp_dip_size) {
  message("BWPS without state-space reduction")
  N_smp <- nrow(obs_smp)
  L_smp <- ncol(obs_smp)
  
  dip <- array(NA, dim = c(N_smp, 2, L_smp))
  for (i in 1:N_smp) {
    print(paste("individual:", i))
    fwd_prob <- cmprunFwdProcedure_Diplo_Std_Pnl(obs_smp[i, ], ref_pnl, theta, lambda)
    smp_dip <- cmprunBkdProcedure_Diplo_Std_Pnl(obs_smp[i, ], ref_pnl, theta, lambda, fwd_prob, smp_dip_size)
    smp_dip_dist <- cmpconstructHist_Diplo(smp_dip, smp_dip_size)
    dip[i, , ] <- smp_dip_dist$diplotype[which(smp_dip_dist$count == max(smp_dip_dist$count))[1], , ]
    
    save(dip, file = "BWPS_Diplo_Std_Pnl_N1000_L0500.rda")
  }
  
  return(dip)
}
#' Compiled version
cmprunBWPS_Diplo_Std_Pnl <- cmpfun(runBWPS_Diplo_Std_Pnl)

##################################################

#' Haplotype phasing with state-space reduction

#' Preprocessing

#' Partition the haplotype reference panel into several genomic blocks with a predetermined block length
#' @param ref_pnl the reference panel
#' @param bck_len the length of the genomic block
#' @param bck_idx the index of the genomic block
#' @return the haplotype reference blocks returned in a list

#' Standard version
partitionRefPanel <- function(ref_pnl, bck_len, bck_idx = NULL) {
  L_pnl <- ncol(ref_pnl)
  
  if (is.null(bck_idx)) {
    bck_idx <- seq(1, L_pnl, bck_len - 1)
    if (bck_idx[length(bck_idx)] < L_pnl) {
      bck_idx <- c(bck_idx, L_pnl)
    }
  }
  
  ref_bck <- list()
  for (i in seq(length(bck_idx) - 1)) {
    ref_bck[[length(ref_bck) + 1]] <- ref_pnl[, bck_idx[i]:bck_idx[i + 1]]
  }
  
  return(lapply(ref_bck, as.matrix))
}
#' Compiled version
cmppartitionRefPanel <- cmpfun(partitionRefPanel)

#########################

#' Collapse the haplotype reference block (i.e., remove the duplicates of the haplotypes in the reference block)
#' @param uncollapsed_ref_bck the uncollapsed haplotype reference block
#' @return the collapsed haplotype reference block with the associated map to the uncollapsed reference block returned in a list

#' Standard version
collapseRefBlock <- function(uncollapsed_ref_bck) {
  unq <- unique(uncollapsed_ref_bck)
  
  env <- new.env()
  env$idx <- list()
  
  cnt <- vector()
  cnt <- apply(unq, 1, 
               function(x) {
                 collapsed_idx <- length(env$idx) + 1
                 uncollapsed_idx <- which(apply((t(t(uncollapsed_ref_bck) - x)) == 0, 1, all))
                 env$idx[[collapsed_idx]] <- uncollapsed_idx
                 env$map[[collapsed_idx]] <- data.frame(collapsed = rep(collapsed_idx, length(uncollapsed_idx)), uncollapsed = uncollapsed_idx)
                 return(length(uncollapsed_idx))
               })
  
  map <- rbindlist(env$map)
  setkey(map, "uncollapsed")
  
  return(list(unq = unq, cnt = cnt, idx = env$idx, map = map))
}
#' Complied version
cmpcollapseRefBlock <- cmpfun(collapseRefBlock)

#########################

#' Partition the switch rates across the chromosome into several genomic blocks with a predetermined block length
#' @param theta the switch rates across the chromosome
#' @param bck_len the length of the genomic block
#' @param bck_idx the index of the genomic block
#' @return the switch rate blocks returned in a list

#' Standard version
partitionSwitchRate <- function(theta, bck_len, bck_idx = NULL) {
  L_pnl <- length(theta) + 1
  
  if (is.null(bck_idx)) {
    bck_idx <- seq(1, L_pnl, bck_len - 1)
    if (bck_idx[length(bck_idx)] < L_pnl) {
      bck_idx <- c(bck_idx, L_pnl)
    }
  }
  
  theta_bck <- list()
  for (i in seq(length(bck_idx) - 1)) {
    theta_bck[[length(theta_bck) + 1]] <- theta[bck_idx[i]:(bck_idx[i + 1] - 1)]
  }
  
  return(lapply(theta_bck, as.vector))
}
#' Compiled version
cmppartitionSwitchRate <- cmpfun(partitionSwitchRate)

#########################

#' Partition the observed genotype across the chromosome into several genomic blocks with a predetermined block length
#' @param obs_gen the observed genotype across the chromosome
#' @param bck_len the length of the genomic block
#' @param bck_idx the index of the genomic block
#' @return the observed genotype blocks returned in a list

#' Standard version
partitionObsGenotype <- function(obs_gen, bck_len, bck_idx = NULL) {
  L_hap <- length(obs_gen)
  
  if (is.null(bck_idx)) {
    bck_idx <- seq(1, L_hap, bck_len - 1)
    if (bck_idx[length(bck_idx)] < L_hap) {
      bck_idx <- c(bck_idx, L_hap)
    }
  }
  
  obs_bck <- list()
  for (i in seq(length(bck_idx) - 1)) {
    obs_bck[[length(obs_bck) + 1]] <- obs_gen[bck_idx[i]:bck_idx[i + 1]]
  }
  
  return(lapply(obs_bck, as.vector))
}
#' Compiled version
cmppartitionObsGenotype <- cmpfun(partitionObsGenotype)

##################################################

#' Calculate the transition probabilities at a specific marker of the chromosome (for collapsed state space) using Eq.20
#' @param ref_hap_a_idx the index of the reference haplotype a
#' @param ref_hap_b_idx the index of the reference haplotype b
#' @param theta the switch rate
#' @param n_bck the number of the duplicates of the each haplotype in the genomic block
#' @param N_bck the number of the haplotypes in the genomic block
#' @param N_pnl the number of the haplotypes in the reference panel
#' @return the transition probabilities at a specific marker of the chromosome returned in a matrix (row = haplotype a, col = haplotype b)

#' Standard version
calculateTransProb_Diplo_Red_Pnl <- function(ref_hap_a_idx, ref_hap_b_idx, theta, n_bck, N_bck, N_pnl) {
  trans_prob <- matrix((theta * n_bck[ref_hap_a_idx] / N_pnl) * (theta * n_bck[ref_hap_b_idx] / N_pnl), nrow = N_bck, ncol = N_bck)
  trans_prob[ref_hap_a_idx, ] <- (theta * n_bck[ref_hap_a_idx] / N_pnl + (1 - theta)) * (theta * n_bck[ref_hap_b_idx] / N_pnl)
  trans_prob[, ref_hap_b_idx] <- (theta * n_bck[ref_hap_a_idx] / N_pnl) * (theta * n_bck[ref_hap_b_idx] / N_pnl + (1 - theta))
  trans_prob[ref_hap_a_idx, ref_hap_b_idx] <- (theta * n_bck[ref_hap_a_idx] / N_pnl + (1 - theta)) * (theta * n_bck[ref_hap_b_idx] / N_pnl + (1 - theta))
  
  return(trans_prob)
}
#' Compiled version
cmpcalculateTransProb_Diplo_Red_Pnl <- cmpfun(calculateTransProb_Diplo_Red_Pnl)

#########################

#' Run the forward procedure with state-space reduction using Eqs.13-59 (no approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param K_bck the number of the blocks
#' @return the collapsed forward probabilities across each genomic block including the uncollapsed forward probabilities at the boundary of each genomic block returned in a list

#' Standard version
runFwdProcedure_Diplo_Red_Pnl_NA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck) {
  message("forward procedure")
  # run the forward procedure with state-space reduction for the genomic block
  runFwdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, N_pnl, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
    collapseFwdProb_Diplo_FullRed_Bck <- function(uncollapsed_fwd_prob_lbndry, collapsed_info) {
      N_bck <- length(collapsed_info$cnt)
      
      collapsed_fwd_prob_lbndry <- matrix(NA, nrow = N_bck, ncol = N_bck)
      for (u in seq(N_bck)) {
        for (v in seq(N_bck)) {
          collapsed_fwd_prob_lbndry[u, v] <- sum(uncollapsed_fwd_prob_lbndry[collapsed_info$idx[[u]], collapsed_info$idx[[v]]])
        }
      }
      
      return(collapsed_fwd_prob_lbndry)
    }
    cmpcollapseFwdProb_Diplo_FullRed_Bck <- cmpfun(collapseFwdProb_Diplo_FullRed_Bck)
    
    updateFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_pre, emn_prob, theta, n_bck, N_pnl, mkr) {
      if (mkr == 2) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre
        }
        collapsed_fwd_prob_nr_re <- function() { 
          apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                function(x) x * rowSums(collapsed_fwd_prob_pre))
        }
        collapsed_fwd_prob_re_nr <- function() {
          apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
      }
      if (mkr >= 3) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
        }
        collapsed_fwd_prob_nr_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_re))
        }
        
        collapsed_fwd_prob_re_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_nr + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$re_nr)), 2, 
                  function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$re_nr + collapsed_fwd_prob_pre$re_re)) + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_re + collapsed_fwd_prob_pre$re_re)), 2, 
                  function(x) x * n_bck) + 
            apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre$total)), 2, 
                  function(x) x * n_bck)
        }
      }
      collapsed_fwd_prob <- list()
      collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
      collapsed_fwd_prob[["nr_re"]] <- collapsed_fwd_prob_nr_re()
      collapsed_fwd_prob[["re_nr"]] <- collapsed_fwd_prob_re_nr()
      collapsed_fwd_prob[["re_re"]] <- collapsed_fwd_prob_re_re()
      collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
      
      return(collapsed_fwd_prob)
    }
    cmpupdateFwdProb_Diplo_FullRed_Bck <- cmpfun(updateFwdProb_Diplo_FullRed_Bck)
    
    uncollapsedFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_nr_nr_rbndry, collapsed_fwd_prob_re_re_rbndry, 
                                                     collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) {
      uncollapsed_fwd_prob_rbndry <- matrix(NA, nrow = N_pnl, ncol = N_pnl)
      for (i in seq(N_pnl)) {
        for (j in seq(N_pnl)) {
          u <- collapsed_info$map[uncollapsed == i]$collapsed
          cnt_bck_hap_a <- collapsed_info$cnt[u]
          v <- collapsed_info$map[uncollapsed == j]$collapsed
          cnt_bck_hap_b <- collapsed_info$cnt[v]
          uncollapsed_fwd_prob_rbndry[i, j] <- 
            collapsed_fwd_prob_nr_nr_rbndry[u, v] * uncollapsed_fwd_prob_lbndry[i, j] / collapsed_fwd_prob_lbndry[u, v] + 
            collapsed_fwd_prob_re_re_rbndry[u, v] / cnt_bck_hap_a / cnt_bck_hap_b
        }
      }
      
      return(uncollapsed_fwd_prob_rbndry)
    }
    cmpuncollapsedFwdProb_Diplo_FullRed_Bck <- cmpfun(uncollapsedFwdProb_Diplo_FullRed_Bck)
    
    collapseFwdProb_Diplo_SemiRed_Bck <- function(uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) {
      N_bck <- length(collapsed_info$cnt)
      
      collapsed_fwd_prob_lbndry <- matrix(NA, N_pnl, N_bck)
      for (i in seq(N_pnl)) {
        for (v in seq(N_bck)) {
          collapsed_fwd_prob_lbndry[i, v] <- sum(uncollapsed_fwd_prob_lbndry[i, collapsed_info$idx[[v]]])
        }
      }
      
      return(collapsed_fwd_prob_lbndry)
    }
    cmpcollapseFwdProb_Diplo_SemiRed_Bck <- cmpfun(collapseFwdProb_Diplo_SemiRed_Bck)
    
    updateFwdProb_Diplo_SemiRed_Bck <- function(collapsed_fwd_prob_pre, emn_prob, theta, n_bck, N_pnl, mkr) {
      if (mkr == 2) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre
        }
        collapsed_fwd_prob_nr_re <- function() {
          apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                function(x) x * rowSums(collapsed_fwd_prob_pre))
        }
      }
      if (mkr >= 3) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
        }
        collapsed_fwd_prob_nr_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_re))
        }
      }
      collapsed_fwd_prob <- list()
      collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
      collapsed_fwd_prob[["nr_re"]] <- collapsed_fwd_prob_nr_re()
      
      return(collapsed_fwd_prob)
    }
    cmpupdateFwdProb_Diplo_SemiRed_Bck <- cmpfun(updateFwdProb_Diplo_SemiRed_Bck)
    
    uncollapsedFwdProb_Diplo_SemiRed_Bck <- function(collapsed_fwd_prob_nr_re_rbndry, N_pnl, collapsed_info) {
      uncollapsed_fwd_prob_nr_re_rbndry <- matrix(NA, nrow = N_pnl, ncol = N_pnl)
      for (i in seq(N_pnl)) {
        for (j in seq(N_pnl)) {
          v <- collapsed_info$map[uncollapsed == j]$collapsed
          cnt_bck_hap_b <- collapsed_info$cnt[v]
          uncollapsed_fwd_prob_nr_re_rbndry[i, j] <- collapsed_fwd_prob_nr_re_rbndry[i, v] / cnt_bck_hap_b
        }
      }
      
      return(uncollapsed_fwd_prob_nr_re_rbndry + t(uncollapsed_fwd_prob_nr_re_rbndry))
    }
    cmpuncollapsedFwdProb_Diplo_SemiRed_Bck <- cmpfun(uncollapsedFwdProb_Diplo_SemiRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- as.vector(collapsed_info$cnt)
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    collapsed_map <- as.vector(unlist(collapsed_info$map[, 1]))
    
    collapsed_fwd_prob <- array(NA, dim = c(N_bck, N_bck, L_bck))
    # print(paste("marker:", 1))
    # full_collapsed_fwd_prob <- cmpcollapseFwdProb_Diplo_FullRed_Bck(uncollapsed_fwd_prob_lbndry, collapsed_info)
    full_collapsed_fwd_prob <- collapseFwdProb_Diplo_FullRed_Bck_NA_cpp(uncollapsed_fwd_prob_lbndry, N_pnl, N_bck, collapsed_map)
    # semi_collapsed_fwd_prob <- cmpcollapseFwdProb_Diplo_SemiRed_Bck(uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info)
    semi_collapsed_fwd_prob <- collapseFwdProb_Diplo_SemiRed_Bck_NA_cpp(uncollapsed_fwd_prob_lbndry, N_pnl, N_bck, collapsed_map)    
    collapsed_fwd_prob[, , 1] <- full_collapsed_fwd_prob
    # print(collapsed_fwd_prob[, , 1])
    if (L_bck > 1) {
      for (l in 2:L_bck) {
        # print(paste("marker:", l))
        # full_collapsed_emn_prob <- cmpcalculateEmnProb_Diplo(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_emn_prob <- calculateEmnProb_Diplo_cpp(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_fwd_prob <- cmpupdateFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob, full_collapsed_emn_prob, theta_bck[l - 1], n_bck, N_pnl, l)
        # semi_collapsed_emn_prob <- cmpcalculateEmnProb_Diplo(obs_bck[l], uncollapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        semi_collapsed_emn_prob <- calculateEmnProb_Diplo_cpp(obs_bck[l], uncollapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        semi_collapsed_fwd_prob <- cmpupdateFwdProb_Diplo_SemiRed_Bck(semi_collapsed_fwd_prob, semi_collapsed_emn_prob, theta_bck[l - 1], n_bck, N_pnl, l)
        if (min(full_collapsed_fwd_prob$total, semi_collapsed_fwd_prob$nr_nr, semi_collapsed_fwd_prob$nr_re) < 1e-30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e+20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e+20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e+20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e+20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e+20
          
          semi_collapsed_fwd_prob$nr_nr <- semi_collapsed_fwd_prob$nr_nr * 1e+20
          semi_collapsed_fwd_prob$nr_re <- semi_collapsed_fwd_prob$nr_re * 1e+20
        }
        if (max(full_collapsed_fwd_prob$total, semi_collapsed_fwd_prob$nr_nr, semi_collapsed_fwd_prob$nr_re) > 1e+30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e-20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e-20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e-20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e-20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e-20
          
          semi_collapsed_fwd_prob$nr_nr <- semi_collapsed_fwd_prob$nr_nr * 1e-20
          semi_collapsed_fwd_prob$nr_re <- semi_collapsed_fwd_prob$nr_re * 1e-20
        }
        collapsed_fwd_prob[, , l] <- full_collapsed_fwd_prob$total
        # print(collapsed_fwd_prob[, , l])
      }
      # uncollapsed_fwd_prob_rbndry <- 
      #   cmpuncollapsedFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
      #                                           collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) + 
      #   cmpuncollapsedFwdProb_Diplo_SemiRed_Bck(semi_collapsed_fwd_prob$nr_re, N_pnl, collapsed_info)
      uncollapsed_fwd_prob_rbndry <- 
        uncollapseFwdProb_Diplo_FullRed_Bck_NA_cpp(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
                                                   collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_map, n_bck) + 
        uncollapseFwdProb_Diplo_SemiRed_Bck_NA_cpp(semi_collapsed_fwd_prob$nr_re, N_pnl, collapsed_map, n_bck)
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    } else {
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_lbndry))
    }
  }
  cmprunFwdProcedure_Diplo_Red_Bck <- cmpfun(runFwdProcedure_Diplo_Red_Bck)
  
  emn_prob_tab <- cmpgenerateEmnProbTab(lambda)
  
  # run the forward procedure
  uncollapsed_fwd_prob_bndry <- array(NA, dim = c(N_pnl, N_pnl, K_bck + 1))
  # uncollapsed_fwd_prob_bndry[, , 1] <- cmpcalculateEmnProb_Diplo(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_bndry[, , 1] <- calculateEmnProb_Diplo_cpp(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  collapsed_fwd_prob <- list()
  
  for (k in 1:K_bck) {
    print(paste("block:", k))
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    fwd_procedure_bck <- cmprunFwdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], N_pnl, 
                                                          uncollapsed_fwd_prob_bndry[, , k], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_bck$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[, , k + 1] <- fwd_procedure_bck$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}
#' Compiled version
cmprunFwdProcedure_Diplo_Red_Pnl_NA <- cmpfun(runFwdProcedure_Diplo_Red_Pnl_NA)

#########################

#' Run the backward procedure with state-space reduction using Eqs.60-73 (no approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param L_pnl the number of the markers across the chromosome in the reference panel
#' @param K_bck the number of the blocks
#' @param collapsed_fwd_prob the full collapsed forward probabilities across each genomic block
#' @param uncollapsed_fwd_prob_bndry the uncollapsed forward probabilities at the boundary of each genomic block
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the sampled diplotypes across the chromosome returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBkdProcedure_Diplo_Red_Pnl_NA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                             collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, smp_dip_size) {
  message("backward procedure")
  # run the backward procedure with state-space reduction for the genomic block
  runBkdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, lambda, N_pnl, 
                                            uncollapsed_bkd_path_rbndry, collapsed_fwd_prob, uncollapsed_fwd_prob_lbndry, smp_dip_size) {
    collapseBkdPath_Diplo_FullRed_Bck <- function(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info) {
      collapsed_bkd_path_rbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        collapsed_bkd_path_rbndry[m, 1] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 1]]$collapsed
        collapsed_bkd_path_rbndry[m, 2] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 2]]$collapsed
      }
      
      return(collapsed_bkd_path_rbndry)
    }
    cmpcollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(collapseBkdPath_Diplo_FullRed_Bck)
    
    uncollapseBkdPath_Diplo_FullRed_Bck <- function(collapsed_bkd_path_lbndry, 
                                                    collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info) {
      uncollapsed_bkd_path_lbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        idx_set_hap_a <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 1]]$uncollapsed
        cnt_bck_hap_a <- length(idx_set_hap_a)
        idx_set_hap_b <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 2]]$uncollapsed
        cnt_bck_hap_b <- length(idx_set_hap_b)
        smp_prob <- as.vector(uncollapsed_fwd_prob_lbndry[idx_set_hap_a, idx_set_hap_b] / collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 1], collapsed_bkd_path_lbndry[m, 2]])
        uncollapsed_bkd_path_lbndry_idx <- sample.int(cnt_bck_hap_a * cnt_bck_hap_b, size = 1, replace = TRUE, prob = smp_prob)
        uncollapsed_bkd_path_lbndry[m, 1] <- 
          idx_set_hap_a[if (uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a == 0) cnt_bck_hap_a else uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a]
        uncollapsed_bkd_path_lbndry[m, 2] <- 
          idx_set_hap_b[if (uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a == 0) uncollapsed_bkd_path_lbndry_idx %/% cnt_bck_hap_a else uncollapsed_bkd_path_lbndry_idx %/% cnt_bck_hap_a + 1]
      }
      
      return(uncollapsed_bkd_path_lbndry)
    }
    cmpuncollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(uncollapseBkdPath_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- collapsed_info$cnt
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    
    collapsed_bkd_path <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    
    # print(paste("marker:", L_bck))
    collapsed_bkd_path[, , L_bck] <- cmpcollapseBkdPath_Diplo_FullRed_Bck(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info)
    # smp_dip[, , L_bck] <- matrix(NA, nrow = smp_dip_size, ncol = 2)
    
    for (l in L_bck:2) {
      # print(paste("marker:", l - 1))
      for (m in 1:smp_dip_size) {
        trans_prob <- cmpcalculateTransProb_Diplo_Red_Pnl(collapsed_bkd_path[m, 1, l], collapsed_bkd_path[m, 2, l], theta_bck[l - 1], n_bck, N_bck, N_pnl)
        smp_prob <- as.vector(collapsed_fwd_prob[, , l - 1] * trans_prob / sum(collapsed_fwd_prob[, , l - 1] * trans_prob))
        collapsed_bkd_path_idx <- sample.int(N_bck * N_bck, size = 1, replace = TRUE, prob = smp_prob)
        collapsed_bkd_path[m, , l - 1] <- 
          if (collapsed_bkd_path_idx %% N_bck == 0) c(N_bck, collapsed_bkd_path_idx %/% N_bck) else c(collapsed_bkd_path_idx %% N_bck, collapsed_bkd_path_idx %/% N_bck + 1)
        smp_dip[m, , l - 1] <- cmpgenerateDiplotype_Diplo(obs_bck[l - 1], collapsed_ref_bck[collapsed_bkd_path[m, , l - 1], l - 1], lambda)
      }
      # print(smp_dip[, , l - 1])
    }
    uncollapsed_bkd_path_lbndry <- cmpuncollapseBkdPath_Diplo_FullRed_Bck(collapsed_bkd_path[, , 1], 
                                                                          collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info)
      
    return(list(diplotype = smp_dip[, , 1:(L_bck - 1)], uncollapsed_bkd_path_lbndry = uncollapsed_bkd_path_lbndry))
  }
  cmprunBkdProcedure_Diplo_Red_Bck <- cmpfun(runBkdProcedure_Diplo_Red_Bck)
  
  # run the backward procedure
  uncollapsed_bkd_path_bndry <- array(NA, dim = c(smp_dip_size, 2, K_bck + 1))
  smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  
  smp_prob <- as.vector(uncollapsed_fwd_prob_bndry[, , K_bck + 1] / sum(uncollapsed_fwd_prob_bndry[, , K_bck + 1]))
  uncollapsed_bkd_path_idx <- sample.int(N_pnl * N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  for (m in 1:smp_dip_size) {
    uncollapsed_bkd_path_bndry[m, , K_bck + 1] <- 
      if (uncollapsed_bkd_path_idx[m] %% N_pnl == 0) c(N_pnl, uncollapsed_bkd_path_idx[m] %/% N_pnl) else c(uncollapsed_bkd_path_idx[m] %% N_pnl, uncollapsed_bkd_path_idx[m] %/% N_pnl + 1)
    smp_dip[m, , L_pnl] <- 
      cmpgenerateDiplotype_Diplo(obs_gen_ptn[[K_bck]][ncol(ref_pnl_ptn[[K_bck]])], ref_pnl_ptn[[K_bck]][uncollapsed_bkd_path_bndry[m, , K_bck + 1], ncol(ref_pnl_ptn[[K_bck]])], lambda)
  }
  
  bck_bndry <- L_pnl
  for (k in K_bck:1) {
    print(paste("block:", k))
    L_bck <- ncol(ref_pnl_ptn[[k]])
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    bkd_procedure_bck <- cmprunBkdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], lambda, N_pnl, 
                                                          uncollapsed_bkd_path_bndry[, , k + 1], collapsed_fwd_prob[[k]], uncollapsed_fwd_prob_bndry[, , k], smp_dip_size)
    smp_dip[, , (bck_bndry - L_bck  + 1):(bck_bndry - 1)] <- bkd_procedure_bck$diplotype
    uncollapsed_bkd_path_bndry[, , k] <- bkd_procedure_bck$uncollapsed_bkd_path_lbndry
    bck_bndry <- bck_bndry - L_bck  + 1
  }
  
  return(smp_dip)
}
#' Compiled version
cmprunBkdProcedure_Diplo_Red_Pnl_NA <- cmpfun(runBkdProcedure_Diplo_Red_Pnl_NA)

#########################

#' Run the block-wise path sampling with state-space reduction using Eqs.13-74 (no approximation)
#' @param obs_smp the observed sample
#' @param ref_pnl the reference panel
#' @param bck_len the length of the genomic block
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the estimated diplotypes across the chromosome achieving the maximum a posterior returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBWPS_Diplo_Red_Pnl_NA <- function(obs_smp, ref_pnl, bck_len, theta, lambda, smp_dip_size) {
  message("BWPS with state-space reduction")
  N_smp <- nrow(obs_smp)
  L_smp <- ncol(obs_smp)
  
  # run the preprocessing
  ref_pnl_ptn <- cmppartitionRefPanel(ref_pnl, bck_len, bck_idx = NULL)
  theta_ptn <- cmppartitionSwitchRate(theta, bck_len, bck_idx = NULL)
  
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  K_bck <- length(ref_pnl_ptn)
  
  dip <- array(NA, dim = c(N_smp, 2, L_smp))
  for (i in 1:N_smp) {
    print(paste("individual:", i))
    obs_gen_ptn <- cmppartitionObsGenotype(obs_smp[i, ], bck_len, bck_idx = NULL)
    fwd_prob <- cmprunFwdProcedure_Diplo_Red_Pnl_NA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck)
    smp_dip <- cmprunBkdProcedure_Diplo_Red_Pnl_NA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                                   fwd_prob$collapsed_fwd_prob, fwd_prob$uncollapsed_fwd_prob_bndry, smp_dip_size)
    smp_dip_dist <- cmpconstructHist_Diplo(smp_dip, smp_dip_size)
    dip[i, , ] <- smp_dip_dist$diplotype[which(smp_dip_dist$count == max(smp_dip_dist$count))[1], , ]
    
    save(dip, file = "BWPS_Diplo_Red_Pnl_NA_N1000_L0500.rda")
  }
  return(dip)
}
#' Compiled version
cmprunBWPS_Diplo_Red_Pnl_NA <- cmpfun(runBWPS_Diplo_Red_Pnl_NA)

##################################################

#' Run the forward procedure with state-space reduction using Eqs.13-29 and 79-80 (the first type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param K_bck the number of the blocks
#' @return the collapsed forward probabilities across each genomic block including the uncollapsed forward probabilities at the boundary of each genomic block returned in a list

#' Standard version
runFwdProcedure_Diplo_Red_Pnl_FA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck) {
  message("forward procedure")
  # run the forward procedure with state-space reduction for the genomic block
  runFwdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, N_pnl, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
    collapseFwdProb_Diplo_FullRed_Bck <- function(uncollapsed_fwd_prob_lbndry, collapsed_info) {
      N_bck <- length(collapsed_info$cnt)
      
      collapsed_fwd_prob_lbndry <- matrix(NA, nrow = N_bck, ncol = N_bck)
      for (u in seq(N_bck)) {
        for (v in seq(N_bck)) {
          collapsed_fwd_prob_lbndry[u, v] <- sum(uncollapsed_fwd_prob_lbndry[collapsed_info$idx[[u]], collapsed_info$idx[[v]]])
        }
      }
      
      return(collapsed_fwd_prob_lbndry)
    }
    cmpcollapseFwdProb_Diplo_FullRed_Bck <- cmpfun(collapseFwdProb_Diplo_FullRed_Bck)
    
    updateFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_pre, emn_prob, theta, n_bck, N_pnl, mkr) {
      if (mkr == 2) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre
        }
        collapsed_fwd_prob_nr_re <- function() { 
          apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                function(x) x * rowSums(collapsed_fwd_prob_pre))
        }
        collapsed_fwd_prob_re_nr <- function() {
          apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
      }
      if (mkr >= 3) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
        }
        collapsed_fwd_prob_nr_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_re))
        }
        
        collapsed_fwd_prob_re_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_nr + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$re_nr)), 2, 
                  function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$re_nr + collapsed_fwd_prob_pre$re_re)) + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_re + collapsed_fwd_prob_pre$re_re)), 2, 
                  function(x) x * n_bck) + 
            apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre$total)), 2, 
                  function(x) x * n_bck)
        }
      }
      collapsed_fwd_prob <- list()
      collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
      collapsed_fwd_prob[["nr_re"]] <- collapsed_fwd_prob_nr_re()
      collapsed_fwd_prob[["re_nr"]] <- collapsed_fwd_prob_re_nr()
      collapsed_fwd_prob[["re_re"]] <- collapsed_fwd_prob_re_re()
      collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
      
      return(collapsed_fwd_prob)
    }
    cmpupdateFwdProb_Diplo_FullRed_Bck <- cmpfun(updateFwdProb_Diplo_FullRed_Bck)
    
    uncollapseFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_nr_nr_rbndry, collapsed_fwd_prob_re_re_rbndry, 
                                                    collapsed_fwd_prob_nr_re_rbndry, collapsed_fwd_prob_re_nr_rbndry, 
                                                    collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) {
      marginal_collapsed_fwd_prob_lbndry <- rowSums(collapsed_fwd_prob_lbndry)      
      marginal_uncollapsed_fwd_prob_lbndry <- rowSums(uncollapsed_fwd_prob_lbndry)
      
      uncollapsed_fwd_prob_rbndry <- matrix(NA, nrow = N_pnl, ncol = N_pnl)
      for (i in seq(N_pnl)) {
        for (j in seq(N_pnl)) {
          u <- collapsed_info$map[uncollapsed == i]$collapsed
          cnt_bck_hap_a <- collapsed_info$cnt[u]
          v <- collapsed_info$map[uncollapsed == j]$collapsed
          cnt_bck_hap_b <- collapsed_info$cnt[v]
          uncollapsed_fwd_prob_rbndry[i, j] <- 
            collapsed_fwd_prob_nr_nr_rbndry[u, v] * uncollapsed_fwd_prob_lbndry[i, j] / collapsed_fwd_prob_lbndry[u, v] + 
            collapsed_fwd_prob_nr_re_rbndry[u, v] * (marginal_uncollapsed_fwd_prob_lbndry[i] / marginal_collapsed_fwd_prob_lbndry[u]) / cnt_bck_hap_b + 
            collapsed_fwd_prob_re_nr_rbndry[u, v] / cnt_bck_hap_a * (marginal_uncollapsed_fwd_prob_lbndry[j] / marginal_collapsed_fwd_prob_lbndry[v]) + 
            collapsed_fwd_prob_re_re_rbndry[u, v] / cnt_bck_hap_a / cnt_bck_hap_b
        }
      }
      
      return(uncollapsed_fwd_prob_rbndry)
    }
    cmpuncollapsedFwdProb_Diplo_FullRed_Bck <- cmpfun(uncollapseFwdProb_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- as.vector(collapsed_info$cnt)
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    collapsed_map <- as.vector(unlist(collapsed_info$map[, 1]))
    
    collapsed_fwd_prob <- array(NA, dim = c(N_bck, N_bck, L_bck))
    # print(paste("marker:", 1))
    # full_collapsed_fwd_prob <- cmpcollapseFwdProb_Diplo_FullRed_Bck(uncollapsed_fwd_prob_lbndry, collapsed_info)
    full_collapsed_fwd_prob <- collapseFwdProb_Diplo_FullRed_Bck_FA_cpp(uncollapsed_fwd_prob_lbndry, N_pnl, N_bck, collapsed_map)
    collapsed_fwd_prob[, , 1] <- full_collapsed_fwd_prob
    # print(collapsed_fwd_prob[, , 1])
    if (L_bck > 1) {
      for (l in 2:L_bck) {
        # print(paste("marker:", l))
        # full_collapsed_emn_prob <- cmpcalculateEmnProb_Diplo(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_emn_prob <- calculateEmnProb_Diplo_cpp(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_fwd_prob <- cmpupdateFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob, full_collapsed_emn_prob, theta_bck[l - 1], n_bck, N_pnl, l)
        if (min(full_collapsed_fwd_prob$total) < 1e-30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e+20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e+20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e+20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e+20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e+20
        }
        if (max(full_collapsed_fwd_prob$total) > 1e+30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e-20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e-20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e-20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e-20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e-20
        }
        collapsed_fwd_prob[, , l] <- full_collapsed_fwd_prob$total
        # print(collapsed_fwd_prob[, , l])
      }
      # uncollapsed_fwd_prob_rbndry <- cmpuncollapsedFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
      #                                                                        full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
      #                                                                        collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info)
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Bck_FA_cpp(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
                                                                                full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
                                                                                collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_map, n_bck)
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    } else {
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_lbndry))
    }
  }
  cmprunFwdProcedure_Diplo_Red_Bck <- cmpfun(runFwdProcedure_Diplo_Red_Bck)
  
  emn_prob_tab <- cmpgenerateEmnProbTab(lambda)
  
  # run the forward procedure
  uncollapsed_fwd_prob_bndry <- array(NA, dim = c(N_pnl, N_pnl, K_bck + 1))
  # uncollapsed_fwd_prob_bndry[, , 1] <- cmpcalculateEmnProb_Diplo(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_bndry[, , 1] <- calculateEmnProb_Diplo_cpp(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  collapsed_fwd_prob <- list()
  
  for (k in 1:K_bck) {
    print(paste("block:", k))
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    fwd_procedure_bck <- cmprunFwdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], N_pnl, uncollapsed_fwd_prob_bndry[, , k], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_bck$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[, , k + 1] <- fwd_procedure_bck$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}
#' Compiled version
cmprunFwdProcedure_Diplo_Red_Pnl_FA <- cmpfun(runFwdProcedure_Diplo_Red_Pnl_FA)

#########################

#' Run the backward procedure with state-space reduction using Eqs.60-73 (the first type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param L_pnl the number of the markers across the chromosome in the reference panel
#' @param K_bck the number of the blocks
#' @param collapsed_fwd_prob the full collapsed forward probabilities across each genomic block
#' @param uncollapsed_fwd_prob_bndry the uncollapsed forward probabilities at the boundary of each genomic block
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the sampled diplotypes across the chromosome returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBkdProcedure_Diplo_Red_Pnl_FA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                             collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, smp_dip_size) {
  message("backward procedure")
  # run the backward procedure with state-space reduction for the genomic block
  runBkdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, lambda, N_pnl, 
                                            uncollapsed_bkd_path_rbndry, collapsed_fwd_prob, uncollapsed_fwd_prob_lbndry, smp_dip_size) {
    collapseBkdPath_Diplo_FullRed_Bck <- function(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info) {
      collapsed_bkd_path_rbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        collapsed_bkd_path_rbndry[m, 1] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 1]]$collapsed
        collapsed_bkd_path_rbndry[m, 2] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 2]]$collapsed
      }
      
      return(collapsed_bkd_path_rbndry)
    }
    cmpcollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(collapseBkdPath_Diplo_FullRed_Bck)
    
    uncollapseBkdPath_Diplo_FullRed_Bck <- function(collapsed_bkd_path_lbndry, 
                                                    collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info) {
      uncollapsed_bkd_path_lbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        idx_set_hap_a <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 1]]$uncollapsed
        cnt_bck_hap_a <- length(idx_set_hap_a)
        idx_set_hap_b <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 2]]$uncollapsed
        cnt_bck_hap_b <- length(idx_set_hap_b)
        smp_prob <- as.vector(uncollapsed_fwd_prob_lbndry[idx_set_hap_a, idx_set_hap_b] / collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 1], collapsed_bkd_path_lbndry[m, 2]])
        uncollapsed_bkd_path_lbndry_idx <- sample.int(cnt_bck_hap_a * cnt_bck_hap_b, size = 1, replace = TRUE, prob = smp_prob)
        uncollapsed_bkd_path_lbndry[m, 1] <- 
          idx_set_hap_a[if (uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a == 0) cnt_bck_hap_a else uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a]
        uncollapsed_bkd_path_lbndry[m, 2] <- 
          idx_set_hap_b[if (uncollapsed_bkd_path_lbndry_idx %% cnt_bck_hap_a == 0) uncollapsed_bkd_path_lbndry_idx %/% cnt_bck_hap_a else uncollapsed_bkd_path_lbndry_idx %/% cnt_bck_hap_a + 1]
      }
      
      return(uncollapsed_bkd_path_lbndry)
    }
    cmpuncollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(uncollapseBkdPath_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- collapsed_info$cnt
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    
    collapsed_bkd_path <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    
    # print(paste("marker:", L_bck))
    collapsed_bkd_path[, , L_bck] <- cmpcollapseBkdPath_Diplo_FullRed_Bck(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info)
    # smp_dip[, , L_bck] <- matrix(NA, nrow = smp_dip_size, ncol = 2)
    
    for (l in L_bck:2) {
      # print(paste("marker:", l - 1))
      for (m in 1:smp_dip_size) {
        trans_prob <- cmpcalculateTransProb_Diplo_Red_Pnl(collapsed_bkd_path[m, 1, l], collapsed_bkd_path[m, 2, l], theta_bck[l - 1], n_bck, N_bck, N_pnl)
        smp_prob <- as.vector(collapsed_fwd_prob[, , l - 1] * trans_prob / sum(collapsed_fwd_prob[, , l - 1] * trans_prob))
        collapsed_bkd_path_idx <- sample.int(N_bck * N_bck, size = 1, replace = TRUE, prob = smp_prob)
        collapsed_bkd_path[m, , l - 1] <- 
          if (collapsed_bkd_path_idx %% N_bck == 0) c(N_bck, collapsed_bkd_path_idx %/% N_bck) else c(collapsed_bkd_path_idx %% N_bck, collapsed_bkd_path_idx %/% N_bck + 1)
        smp_dip[m, , l - 1] <- cmpgenerateDiplotype_Diplo(obs_bck[l - 1], collapsed_ref_bck[collapsed_bkd_path[m, , l - 1], l - 1], lambda)
      }
      # print(smp_dip[, , l - 1])
    }
    uncollapsed_bkd_path_lbndry <- cmpuncollapseBkdPath_Diplo_FullRed_Bck(collapsed_bkd_path[, , 1], 
                                                                          collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info)
    
    return(list(diplotype = smp_dip[, , 1:(L_bck - 1)], uncollapsed_bkd_path_lbndry = uncollapsed_bkd_path_lbndry))
  }
  cmprunBkdProcedure_Diplo_Red_Bck <- cmpfun(runBkdProcedure_Diplo_Red_Bck)
  
  # run the backward procedure
  uncollapsed_bkd_path_bndry <- array(NA, dim = c(smp_dip_size, 2, K_bck + 1))
  smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  
  smp_prob <- as.vector(uncollapsed_fwd_prob_bndry[, , K_bck + 1] / sum(uncollapsed_fwd_prob_bndry[, , K_bck + 1]))
  uncollapsed_bkd_path_idx <- sample.int(N_pnl * N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  for (m in 1:smp_dip_size) {
    uncollapsed_bkd_path_bndry[m, , K_bck + 1] <- 
      if (uncollapsed_bkd_path_idx[m] %% N_pnl == 0) c(N_pnl, uncollapsed_bkd_path_idx[m] %/% N_pnl) else c(uncollapsed_bkd_path_idx[m] %% N_pnl, uncollapsed_bkd_path_idx[m] %/% N_pnl + 1)
    smp_dip[m, , L_pnl] <- 
      cmpgenerateDiplotype_Diplo(obs_gen_ptn[[K_bck]][ncol(ref_pnl_ptn[[K_bck]])], ref_pnl_ptn[[K_bck]][uncollapsed_bkd_path_bndry[m, , K_bck + 1], ncol(ref_pnl_ptn[[K_bck]])], lambda)
  }
  
  bck_bndry <- L_pnl
  for (k in K_bck:1) {
    print(paste("block:", k))
    L_bck <- ncol(ref_pnl_ptn[[k]])
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    bkd_procedure_bck <- cmprunBkdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], lambda, N_pnl, 
                                                          uncollapsed_bkd_path_bndry[, , k + 1], collapsed_fwd_prob[[k]], uncollapsed_fwd_prob_bndry[, , k], smp_dip_size)
    smp_dip[, , (bck_bndry - L_bck  + 1):(bck_bndry - 1)] <- bkd_procedure_bck$diplotype
    uncollapsed_bkd_path_bndry[, , k] <- bkd_procedure_bck$uncollapsed_bkd_path_lbndry
    bck_bndry <- bck_bndry - L_bck  + 1
  }
  
  return(smp_dip)
}
#' Compiled version
cmprunBkdProcedure_Diplo_Red_Pnl_FA <- cmpfun(runBkdProcedure_Diplo_Red_Pnl_FA)

#########################

#' Run the block-wise path sampling with state-space reduction using Eqs.13-29, 60-74 and 79-80 (the first type of the approximation)
#' @param obs_smp the observed sample
#' @param ref_pnl the reference panel
#' @param bck_len the length of the genomic block
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the estimated diplotypes across the chromosome achieving the maximum a posterior returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBWPS_Diplo_Red_Pnl_FA <- function(obs_smp, ref_pnl, bck_len, theta, lambda, smp_dip_size) {
  message("BWPS with state-space reduction using the 1st type of the approximation at the boundary of each block")
  N_smp <- nrow(obs_smp)
  L_smp <- ncol(obs_smp)
  
  # run the preprocessing
  ref_pnl_ptn <- cmppartitionRefPanel(ref_pnl, bck_len, bck_idx = NULL)
  theta_ptn <- cmppartitionSwitchRate(theta, bck_len, bck_idx = NULL)
  
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  K_bck <- length(ref_pnl_ptn)
  
  dip <- array(NA, dim = c(N_smp, 2, L_smp))
  for (i in 1:N_smp) {
    print(paste("individual:", i))
    obs_gen_ptn <- cmppartitionObsGenotype(obs_smp[i, ], bck_len, bck_idx = NULL)
    fwd_prob <- cmprunFwdProcedure_Diplo_Red_Pnl_FA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck)
    smp_dip <- cmprunBkdProcedure_Diplo_Red_Pnl_FA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                                   fwd_prob$collapsed_fwd_prob, fwd_prob$uncollapsed_fwd_prob_bndry, smp_dip_size)
    smp_dip_dist <- cmpconstructHist_Diplo(smp_dip, smp_dip_size)
    dip[i, , ] <- smp_dip_dist$diplotype[which(smp_dip_dist$count == max(smp_dip_dist$count))[1], , ]
    
    save(dip, file = "BWPS_Diplo_Red_Pnl_FA_N1000_L0500.rda")
  }
  return(dip)
}
#' Compiled version
cmprunBWPS_Diplo_Red_Pnl_FA <- cmpfun(runBWPS_Diplo_Red_Pnl_FA)

##################################################

#' Run the forward procedure with state-space reduction using Eqs.13-29 and 83-84 (the second type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param K_bck the number of the blocks
#' @return the collapsed forward probabilities across each genomic block including the (marginal) uncollapsed forward probabilities at the boundary of each genomic block returned in a list

#' Standard version
runFwdProcedure_Diplo_Red_Pnl_SA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck) {
  message("forward procedure")
  # run the forward procedure with state-space reduction for the genomic block
  runFwdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, N_pnl, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
    collapseFwdProb_Diplo_FullRed_Bck <- function(uncollapsed_fwd_prob_lbndry, collapsed_info) {
      N_bck <- length(collapsed_info$cnt)
      
      collapsed_fwd_prob_lbndry <- matrix(NA, nrow = N_bck, ncol = N_bck)
      for (u in seq(N_bck)) {
        for (v in seq(N_bck)) {
          collapsed_fwd_prob_lbndry[u, v] <- sum(uncollapsed_fwd_prob_lbndry[collapsed_info$idx[[u]], collapsed_info$idx[[v]]])
        }
      }
      
      return(collapsed_fwd_prob_lbndry)
    }
    cmpcollapseFwdProb_Diplo_FullRed_Bck <- cmpfun(collapseFwdProb_Diplo_FullRed_Bck)
    
    updateFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_pre, emn_prob, theta, n_bck, N_pnl, mkr) {
      if (mkr == 2) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre
        }
        collapsed_fwd_prob_nr_re <- function() { 
          apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                function(x) x * rowSums(collapsed_fwd_prob_pre))
        }
        collapsed_fwd_prob_re_nr <- function() {
          apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
      }
      if (mkr >= 3) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
        }
        collapsed_fwd_prob_nr_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_re))
        }
        
        collapsed_fwd_prob_re_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_nr + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$re_nr)), 2, 
                  function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$re_nr + collapsed_fwd_prob_pre$re_re)) + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_re + collapsed_fwd_prob_pre$re_re)), 2, 
                  function(x) x * n_bck) + 
            apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre$total)), 2, 
                  function(x) x * n_bck)
        }
      }
      collapsed_fwd_prob <- list()
      collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
      collapsed_fwd_prob[["nr_re"]] <- collapsed_fwd_prob_nr_re()
      collapsed_fwd_prob[["re_nr"]] <- collapsed_fwd_prob_re_nr()
      collapsed_fwd_prob[["re_re"]] <- collapsed_fwd_prob_re_re()
      collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
      
      return(collapsed_fwd_prob)
    }
    cmpupdateFwdProb_Diplo_FullRed_Bck <- cmpfun(updateFwdProb_Diplo_FullRed_Bck)
    
    uncollapsedFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_nr_nr_rbndry, collapsed_fwd_prob_re_re_rbndry, 
                                                     collapsed_fwd_prob_nr_re_rbndry, collapsed_fwd_prob_re_nr_rbndry, 
                                                     collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) {
      marginal_collapsed_fwd_prob_lbndry <- rowSums(collapsed_fwd_prob_lbndry)
      marginal_uncollapsed_fwd_prob_lbndry <- rowSums(uncollapsed_fwd_prob_lbndry)
      
      uncollapsed_fwd_prob_rbndry <- matrix(NA, nrow = N_pnl, ncol = N_pnl)
      for (i in seq(N_pnl)) {
        for (j in seq(N_pnl)) {
          u <- collapsed_info$map[uncollapsed == i]$collapsed
          cnt_bck_hap_a <- collapsed_info$cnt[u]
          v <- collapsed_info$map[uncollapsed == j]$collapsed
          cnt_bck_hap_b <- collapsed_info$cnt[v]
          uncollapsed_fwd_prob_rbndry[i, j] <- 
            collapsed_fwd_prob_nr_nr_rbndry[u, v] * (marginal_uncollapsed_fwd_prob_lbndry[i] / marginal_collapsed_fwd_prob_lbndry[u]) * (marginal_uncollapsed_fwd_prob_lbndry[j] / marginal_collapsed_fwd_prob_lbndry[v]) + 
            collapsed_fwd_prob_nr_re_rbndry[u, v] * (marginal_uncollapsed_fwd_prob_lbndry[i] / marginal_collapsed_fwd_prob_lbndry[u]) / cnt_bck_hap_b + 
            collapsed_fwd_prob_re_nr_rbndry[u, v] / cnt_bck_hap_a * (marginal_uncollapsed_fwd_prob_lbndry[j] / marginal_collapsed_fwd_prob_lbndry[v]) + 
            collapsed_fwd_prob_re_re_rbndry[u, v] / cnt_bck_hap_a / cnt_bck_hap_b
        }
      }
      
      return(uncollapsed_fwd_prob_rbndry)
    }
    cmpuncollapsedFwdProb_Diplo_FullRed_Bck <- cmpfun(uncollapsedFwdProb_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- as.vector(collapsed_info$cnt)
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    collapsed_map <- as.vector(unlist(collapsed_info$map[, 1]))
    
    collapsed_fwd_prob <- array(NA, dim = c(N_bck, N_bck, L_bck))
    # print(paste("marker:", 1))
    # full_collapsed_fwd_prob <- cmpcollapseFwdProb_Diplo_FullRed_Bck(uncollapsed_fwd_prob_lbndry, collapsed_info)
    full_collapsed_fwd_prob <- collapseFwdProb_Diplo_FullRed_Bck_SA_cpp(uncollapsed_fwd_prob_lbndry, N_pnl, N_bck, collapsed_map)
    collapsed_fwd_prob[, , 1] <- full_collapsed_fwd_prob
    # print(collapsed_fwd_prob[, , 1])
    if (L_bck > 1) {
      for (l in 2:L_bck) {
        # print(paste("marker:", l))
        # full_collapsed_emn_prob <- cmpcalculateEmnProb_Diplo(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_emn_prob <- calculateEmnProb_Diplo_cpp(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_fwd_prob <- cmpupdateFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob, full_collapsed_emn_prob, theta_bck[l - 1], n_bck, N_pnl, l)
        if (min(full_collapsed_fwd_prob$total) < 1e-30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e+20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e+20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e+20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e+20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e+20
        }
        if (max(full_collapsed_fwd_prob$total) > 1e+30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e-20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e-20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e-20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e-20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e-20
        }
        collapsed_fwd_prob[, , l] <- full_collapsed_fwd_prob$total
        # print(collapsed_fwd_prob[, , l])
      }
      # uncollapsed_fwd_prob_rbndry <- cmpuncollapsedFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
      #                                                                        full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
      #                                                                        collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info)
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Bck_SA_cpp(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
                                                                                full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
                                                                                collapsed_fwd_prob[, , 1], uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_map, n_bck)
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    } else {
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_lbndry))
    }
  }
  cmprunFwdProcedure_Diplo_Red_Bck <- cmpfun(runFwdProcedure_Diplo_Red_Bck)
  
  emn_prob_tab <- cmpgenerateEmnProbTab(lambda)
  
  # run the forward procedure
  uncollapsed_fwd_prob_bndry <- matrix(NA, nrow = N_pnl, ncol = K_bck + 1)
  # uncollapsed_fwd_prob_lbndry <- cmpcalculateEmnProb_Diplo(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_lbndry <- calculateEmnProb_Diplo_cpp(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_bndry[, 1] <- rowSums(uncollapsed_fwd_prob_lbndry)
  collapsed_fwd_prob <- list()
  
  for (k in seq(K_bck)) {
    print(paste("block:", k))
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    fwd_procedure_bck <- cmprunFwdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], N_pnl, uncollapsed_fwd_prob_lbndry, emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_bck$collapsed_fwd_prob
    uncollapsed_fwd_prob_lbndry <- fwd_procedure_bck$uncollapsed_fwd_prob_rbndry
    uncollapsed_fwd_prob_bndry[, k + 1] <- rowSums(uncollapsed_fwd_prob_lbndry)
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}
#' Compiled version
cmprunFwdProcedure_Diplo_Red_Pnl_SA <- cmpfun(runFwdProcedure_Diplo_Red_Pnl_SA)

#########################

#' Run the backward procedure with state-space reduction using Eqs.60-73 and 85-86 (the second type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param L_pnl the number of the markers across the chromosome in the reference panel
#' @param K_bck the number of the blocks
#' @param collapsed_fwd_prob the full collapsed forward probabilities across each genomic block
#' @param uncollapsed_fwd_prob_bndry the uncollapsed (marginal) forward probabilities at the boundary of each genomic block
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the sampled diplotypes across the chromosome returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBkdProcedure_Diplo_Red_Pnl_SA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                             collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, smp_dip_size) {
  message("backward procedure")
  # run the backward procedure with state-space reduction for the genomic block
  runBkdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, lambda, N_pnl, 
                                            uncollapsed_bkd_path_rbndry, collapsed_fwd_prob, marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size) {
    collapseBkdPath_Diplo_FullRed_Bck <- function(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info) {
      collapsed_bkd_path_rbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        collapsed_bkd_path_rbndry[m, 1] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 1]]$collapsed
        collapsed_bkd_path_rbndry[m, 2] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 2]]$collapsed
      }
      
      return(collapsed_bkd_path_rbndry)
    }
    cmpcollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(collapseBkdPath_Diplo_FullRed_Bck)
    
    uncollapseBkdPath_Diplo_FullRed_Bck <- function(collapsed_bkd_path_lbndry, 
                                                    collapsed_fwd_prob_lbndry, marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info) {
      marginal_collapsed_fwd_prob_lbndry <- rowSums(collapsed_fwd_prob_lbndry)
      
      uncollapsed_bkd_path_lbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        idx_set_hap_a <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 1]]$uncollapsed
        cnt_bck_hap_a <- length(idx_set_hap_a)
        smp_prob <- as.vector(marginal_uncollapsed_fwd_prob_lbndry[idx_set_hap_a] / marginal_collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 1]])
        uncollapsed_bkd_path_lbndry[m, 1] <- idx_set_hap_a[sample.int(cnt_bck_hap_a, size = 1, replace = TRUE, prob = smp_prob)]
        
        idx_set_hap_b <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 2]]$uncollapsed
        cnt_bck_hap_b <- length(idx_set_hap_b)
        smp_prob <- as.vector(marginal_uncollapsed_fwd_prob_lbndry[idx_set_hap_b] / marginal_collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 2]])
        uncollapsed_bkd_path_lbndry[m, 2] <- idx_set_hap_b[sample.int(cnt_bck_hap_b, size = 1, replace = TRUE, prob = smp_prob)]
      }
      
      return(uncollapsed_bkd_path_lbndry)
    }
    cmpuncollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(uncollapseBkdPath_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- collapsed_info$cnt
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    
    collapsed_bkd_path <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    
    # print(paste("marker:", L_bck))
    collapsed_bkd_path[, , L_bck] <- cmpcollapseBkdPath_Diplo_FullRed_Bck(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info)
    # smp_dip[, , L_bck] <- matrix(NA, nrow = smp_dip_size, ncol = 2)
    
    for (l in L_bck:2) {
      # print(paste("marker:", l - 1))
      for (m in 1:smp_dip_size) {
        trans_prob <- cmpcalculateTransProb_Diplo_Red_Pnl(collapsed_bkd_path[m, 1, l], collapsed_bkd_path[m, 2, l], theta_bck[l - 1], n_bck, N_bck, N_pnl)
        smp_prob <- as.vector(collapsed_fwd_prob[, , l - 1] * trans_prob / sum(collapsed_fwd_prob[, , l - 1] * trans_prob))
        collapsed_bkd_path_idx <- sample.int(N_bck * N_bck, size = 1, replace = TRUE, prob = smp_prob)
        collapsed_bkd_path[m, , l - 1] <- 
          if (collapsed_bkd_path_idx %% N_bck == 0) c(N_bck, collapsed_bkd_path_idx %/% N_bck) else c(collapsed_bkd_path_idx %% N_bck, collapsed_bkd_path_idx %/% N_bck + 1)
        smp_dip[m, , l - 1] <- cmpgenerateDiplotype_Diplo(obs_bck[l - 1], collapsed_ref_bck[collapsed_bkd_path[m, , l - 1], l - 1], lambda)
      }
      # print(smp_dip[, , l - 1])
    }
    uncollapsed_bkd_path_lbndry <- cmpuncollapseBkdPath_Diplo_FullRed_Bck(collapsed_bkd_path[, , 1], 
                                                                          collapsed_fwd_prob[, , 1], marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info)
    
    return(list(diplotype = smp_dip[, , 1:(L_bck - 1)], uncollapsed_bkd_path_lbndry = uncollapsed_bkd_path_lbndry))
  }
  cmprunBkdProcedure_Diplo_Red_Bck <- cmpfun(runBkdProcedure_Diplo_Red_Bck)
  
  # run the backward procedure
  uncollapsed_bkd_path_bndry <- array(NA, dim = c(smp_dip_size, 2, K_bck + 1))
  smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  
  smp_prob <- as.vector(uncollapsed_fwd_prob_bndry[, K_bck + 1] / sum(uncollapsed_fwd_prob_bndry[, K_bck + 1]))
  uncollapsed_bkd_path_bndry[, 1, K_bck + 1] <- sample.int(N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  uncollapsed_bkd_path_bndry[, 2, K_bck + 1] <- sample.int(N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  for (m in 1:smp_dip_size) {
    smp_dip[m, , L_pnl] <- 
      cmpgenerateDiplotype_Diplo(obs_gen_ptn[[K_bck]][ncol(ref_pnl_ptn[[K_bck]])], ref_pnl_ptn[[K_bck]][uncollapsed_bkd_path_bndry[m, , K_bck + 1], ncol(ref_pnl_ptn[[K_bck]])], lambda)
  }
  
  bck_bndry <- L_pnl
  for (k in K_bck:1) {
    print(paste("block:", k))
    L_bck <- ncol(ref_pnl_ptn[[k]])
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    bkd_procedure_bck <- cmprunBkdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], lambda, N_pnl, 
                                                          uncollapsed_bkd_path_bndry[, , k + 1], collapsed_fwd_prob[[k]], uncollapsed_fwd_prob_bndry[, k], smp_dip_size)
    smp_dip[, , (bck_bndry - L_bck  + 1):(bck_bndry - 1)] <- bkd_procedure_bck$diplotype
    uncollapsed_bkd_path_bndry[, , k] <- bkd_procedure_bck$uncollapsed_bkd_path_lbndry
    bck_bndry <- bck_bndry - L_bck  + 1
  }
  
  return(smp_dip)
}
#' Compiled version
cmprunBkdProcedure_Diplo_Red_Pnl_SA <- cmpfun(runBkdProcedure_Diplo_Red_Pnl_SA)

#########################

#' Run the block-wise path sampling with state-space reduction using Eqs.13-29, 60-74 and 83-86 (the second type of the approximation)
#' @param obs_smp the observed sample
#' @param ref_pnl the reference panel
#' @param bck_len the length of the genomic block
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the estimated diplotypes across the chromosome achieving the maximum a posterior returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBWPS_Diplo_Red_Pnl_SA <- function(obs_smp, ref_pnl, bck_len, theta, lambda, smp_dip_size) {
  message("BWPS with state-space reduction using the 2nd type of the approximation at the boundary of each block")
  N_smp <- nrow(obs_smp)
  L_smp <- ncol(obs_smp)
  
  # run the preprocessing
  ref_pnl_ptn <- cmppartitionRefPanel(ref_pnl, bck_len, bck_idx = NULL)
  theta_ptn <- cmppartitionSwitchRate(theta, bck_len, bck_idx = NULL)
  
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  K_bck <- length(ref_pnl_ptn)
  
  dip <- array(NA, dim = c(N_smp, 2, L_smp))
  for (i in 1:N_smp) {
    print(paste("individual:", i))
    obs_gen_ptn <- cmppartitionObsGenotype(obs_smp[i, ], bck_len, bck_idx = NULL)
    fwd_prob <- cmprunFwdProcedure_Diplo_Red_Pnl_SA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck)
    smp_dip <- cmprunBkdProcedure_Diplo_Red_Pnl_SA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                                   fwd_prob$collapsed_fwd_prob, fwd_prob$uncollapsed_fwd_prob_bndry, smp_dip_size)
    smp_dip_dist <- cmpconstructHist_Diplo(smp_dip, smp_dip_size)
    dip[i, , ] <- smp_dip_dist$diplotype[which(smp_dip_dist$count == max(smp_dip_dist$count))[1], , ]
    
    save(dip, file = "BWPS_Diplo_Red_Pnl_SA_N1000_L0500.rda")
  }
  return(dip)
}
#' Compiled version
cmprunBWPS_Diplo_Red_Pnl_SA <- cmpfun(runBWPS_Diplo_Red_Pnl_SA)

##################################################

#' Run the forward procedure with state-space reduction using Eqs.13-29 and 89-92 (the third type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param K_bck the number of the blocks
#' @return the collapsed forward probabilities across each genomic block including the (marginal) uncollapsed forward probabilities at the boundary of each genomic block returned in a list

#' Standard version
runFwdProcedure_Diplo_Red_Pnl_TA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck) {
  message("forward procedure")
  # run the forward procedure with state-space reduction for the genomic block
  runFwdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, N_pnl, marginal_uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
    collapseFwdProb_Diplo_FullRed_Bck <- function(marginal_uncollapsed_fwd_prob_lbndry, collapsed_info) {
      N_bck <- length(collapsed_info$cnt)
      
      marginal_collapsed_fwd_prob_lbndry <- rep(NA, N_bck)
      for (u in seq(N_bck)) {
        marginal_collapsed_fwd_prob_lbndry[u] <- sum(marginal_uncollapsed_fwd_prob_lbndry[collapsed_info$idx[[u]]])
      }
      
      collapsed_fwd_prob_lbndry <- marginal_collapsed_fwd_prob_lbndry %*% t(marginal_collapsed_fwd_prob_lbndry)
      
      return(collapsed_fwd_prob_lbndry)
    }
    cmpcollapseFwdProb_Diplo_FullRed_Bck <- cmpfun(collapseFwdProb_Diplo_FullRed_Bck)
    
    updateFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_pre, emn_prob, theta, n_bck, N_pnl, mkr) {
      if (mkr == 2) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre
        }
        collapsed_fwd_prob_nr_re <- function() { 
          apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                function(x) x * rowSums(collapsed_fwd_prob_pre))
        }
        collapsed_fwd_prob_re_nr <- function() {
          apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre)), 2, 
                function(x) x * n_bck)
        }
      }
      if (mkr >= 3) {
        collapsed_fwd_prob_nr_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
        }
        collapsed_fwd_prob_nr_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_re))
        }
        
        collapsed_fwd_prob_re_nr <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_nr + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$re_nr)), 2, 
                  function(x) x * n_bck)
        }
        collapsed_fwd_prob_re_re <- function() {
          emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$re_re + 
            apply(emn_prob %*% diag((1 - theta) * (theta * n_bck / N_pnl)), 2, 
                  function(x) x * rowSums(collapsed_fwd_prob_pre$re_nr + collapsed_fwd_prob_pre$re_re)) + 
            apply(emn_prob %*% diag((theta / N_pnl) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_re + collapsed_fwd_prob_pre$re_re)), 2, 
                  function(x) x * n_bck) + 
            apply(emn_prob %*% diag((theta * n_bck / N_pnl) * (theta / N_pnl) * sum(collapsed_fwd_prob_pre$total)), 2, 
                  function(x) x * n_bck)
        }
      }
      collapsed_fwd_prob <- list()
      collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
      collapsed_fwd_prob[["nr_re"]] <- collapsed_fwd_prob_nr_re()
      collapsed_fwd_prob[["re_nr"]] <- collapsed_fwd_prob_re_nr()
      collapsed_fwd_prob[["re_re"]] <- collapsed_fwd_prob_re_re()
      collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
      
      return(collapsed_fwd_prob)
    }
    cmpupdateFwdProb_Diplo_FullRed_Bck <- cmpfun(updateFwdProb_Diplo_FullRed_Bck)
    
    uncollapsedFwdProb_Diplo_FullRed_Bck <- function(collapsed_fwd_prob_nr_nr_rbndry, collapsed_fwd_prob_re_re_rbndry, 
                                                     collapsed_fwd_prob_nr_re_rbndry, collapsed_fwd_prob_re_nr_rbndry, 
                                                     collapsed_fwd_prob_lbndry, marginal_uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info) {
      marginal_collapsed_fwd_prob_lbndry <- rowSums(collapsed_fwd_prob_lbndry)
      
      marginal_uncollapsed_fwd_prob_rbndry <- rep(NA, N_pnl)
      for (i in seq(N_pnl)) {
        u <- collapsed_info$map[uncollapsed == i]$collapsed
        cnt_bck_hap_a <- collapsed_info$cnt[u]
        marginal_uncollapsed_fwd_prob_rbndry[i] <- 
          sum(collapsed_fwd_prob_nr_nr_rbndry[u, ] + collapsed_fwd_prob_nr_re_rbndry[u, ]) * (marginal_uncollapsed_fwd_prob_lbndry[i] / marginal_collapsed_fwd_prob_lbndry[u]) + 
          sum(collapsed_fwd_prob_re_nr_rbndry[u, ] + collapsed_fwd_prob_re_re_rbndry[u, ]) / cnt_bck_hap_a
      }
      
      return(marginal_uncollapsed_fwd_prob_rbndry)
    }
    cmpuncollapsedFwdProb_Diplo_FullRed_Bck <- cmpfun(uncollapsedFwdProb_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- as.vector(collapsed_info$cnt)
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    collapsed_map <- as.vector(unlist(collapsed_info$map[, 1]))
    
    collapsed_fwd_prob <- array(NA, dim = c(N_bck, N_bck, L_bck))
    # print(paste("marker:", 1))
    # full_collapsed_fwd_prob <- cmpcollapseFwdProb_Diplo_FullRed_Bck(marginal_uncollapsed_fwd_prob_lbndry, collapsed_info)
    full_collapsed_fwd_prob <- collapseFwdProb_Diplo_FullRed_Bck_TA_cpp(marginal_uncollapsed_fwd_prob_lbndry, N_pnl, N_bck, collapsed_map)
    collapsed_fwd_prob[, , 1] <- full_collapsed_fwd_prob
    # print(collapsed_fwd_prob[, , 1])
    if (L_bck > 1) {
      for (l in 2:L_bck) {
        # print(paste("marker:", l))
        # full_collapsed_emn_prob <- cmpcalculateEmnProb_Diplo(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_emn_prob <- calculateEmnProb_Diplo_cpp(obs_bck[l], collapsed_ref_bck[, l], collapsed_ref_bck[, l], emn_prob_tab)
        full_collapsed_fwd_prob <- cmpupdateFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob, full_collapsed_emn_prob, theta_bck[l - 1], n_bck, N_pnl, l)
        if (min(full_collapsed_fwd_prob$total) < 1e-30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e+20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e+20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e+20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e+20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e+20
        }
        if (max(full_collapsed_fwd_prob$total) > 1e+30) {
          full_collapsed_fwd_prob$total <- full_collapsed_fwd_prob$total * 1e-20
          full_collapsed_fwd_prob$nr_nr <- full_collapsed_fwd_prob$nr_nr * 1e-20
          full_collapsed_fwd_prob$nr_re <- full_collapsed_fwd_prob$nr_re * 1e-20
          full_collapsed_fwd_prob$re_nr <- full_collapsed_fwd_prob$re_nr * 1e-20
          full_collapsed_fwd_prob$re_re <- full_collapsed_fwd_prob$re_re * 1e-20
        }
        collapsed_fwd_prob[, , l] <- full_collapsed_fwd_prob$total
        # print(collapsed_fwd_prob[, , l])
      }
      # marginal_uncollapsed_fwd_prob_rbndry <- cmpuncollapsedFwdProb_Diplo_FullRed_Bck(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
      #                                                                                 full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
      #                                                                                 collapsed_fwd_prob[, , 1], marginal_uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_info)
      marginal_uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Bck_TA_cpp(full_collapsed_fwd_prob$nr_nr, full_collapsed_fwd_prob$re_re, 
                                                                                         full_collapsed_fwd_prob$nr_re, full_collapsed_fwd_prob$re_nr, 
                                                                                         collapsed_fwd_prob[, , 1], marginal_uncollapsed_fwd_prob_lbndry, N_pnl, collapsed_map, n_bck)
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = marginal_uncollapsed_fwd_prob_rbndry))
    } else {
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = marginal_uncollapsed_fwd_prob_lbndry))
    }
  }
  cmprunFwdProcedure_Diplo_Red_Bck <- cmpfun(runFwdProcedure_Diplo_Red_Bck)
  
  emn_prob_tab <- cmpgenerateEmnProbTab(lambda)
  
  # run the forward procedure
  uncollapsed_fwd_prob_bndry <- matrix(NA, nrow = N_pnl, ncol = K_bck + 1)
  # uncollapsed_fwd_prob_lbndry <- cmpcalculateEmnProb_Diplo(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_lbndry <- calculateEmnProb_Diplo_cpp(obs_gen_ptn[[1]][1], ref_pnl_ptn[[1]][, 1], ref_pnl_ptn[[1]][, 1], emn_prob_tab) / (N_pnl ^ 2)
  uncollapsed_fwd_prob_bndry[, 1] <- rowSums(uncollapsed_fwd_prob_lbndry)
  collapsed_fwd_prob <- list()
  
  for (k in seq(K_bck)) {
    print(paste("block:", k))
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    fwd_procedure_bck <- cmprunFwdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], N_pnl, uncollapsed_fwd_prob_bndry[, k], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_bck$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[, k + 1] <- fwd_procedure_bck$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}
#' Compiled version
cmprunFwdProcedure_Diplo_Red_Pnl_TA <- cmpfun(runFwdProcedure_Diplo_Red_Pnl_TA)

#########################

#' Run the backward procedure with state-space reduction using Eqs.60-73 and 85-86 (the third type of the approximation)
#' @param obs_gen_ptn the observed genotype blocks
#' @param ref_pnl_ptn the reference panel blocks
#' @param theta_ptn the switch rate blocks
#' @param lambda the mutation rate
#' @param N_pnl the number of the haplotypes in the reference panel
#' @param L_pnl the number of the markers across the chromosome in the reference panel
#' @param K_bck the number of the blocks
#' @param collapsed_fwd_prob the full collapsed forward probabilities across each genomic block
#' @param uncollapsed_fwd_prob_bndry the uncollapsed (marginal) forward probabilities at the boundary of each genomic block
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the sampled diplotypes across the chromosome returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBkdProcedure_Diplo_Red_Pnl_TA <- function(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                             collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, smp_dip_size) {
  message("backward procedure")
  # run the backward procedure with state-space reduction for the genomic block
  runBkdProcedure_Diplo_Red_Bck <- function(obs_bck, uncollapsed_ref_bck, collapsed_info, theta_bck, lambda, N_pnl, 
                                            uncollapsed_bkd_path_rbndry, collapsed_fwd_prob, marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size) {
    collapseBkdPath_Diplo_FullRed_Bck <- function(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info) {
      collapsed_bkd_path_rbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        collapsed_bkd_path_rbndry[m, 1] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 1]]$collapsed
        collapsed_bkd_path_rbndry[m, 2] <- collapsed_info$map[uncollapsed == uncollapsed_bkd_path_rbndry[m, 2]]$collapsed
      }
      
      return(collapsed_bkd_path_rbndry)
    }
    cmpcollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(collapseBkdPath_Diplo_FullRed_Bck)
    
    uncollapseBkdPath_Diplo_FullRed_Bck <- function(collapsed_bkd_path_lbndry, 
                                                    collapsed_fwd_prob_lbndry, marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info) {
      marginal_collapsed_fwd_prob_lbndry <- rowSums(collapsed_fwd_prob_lbndry)

      uncollapsed_bkd_path_lbndry <- matrix(NA, nrow = smp_dip_size, ncol = 2)
      for (m in 1:smp_dip_size) {
        idx_set_hap_a <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 1]]$uncollapsed
        cnt_bck_hap_a <- length(idx_set_hap_a)
        smp_prob <- as.vector(marginal_uncollapsed_fwd_prob_lbndry[idx_set_hap_a] / marginal_collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 1]])
        uncollapsed_bkd_path_lbndry[m, 1] <- idx_set_hap_a[sample.int(cnt_bck_hap_a, size = 1, replace = TRUE, prob = smp_prob)]
        
        idx_set_hap_b <- collapsed_info$map[collapsed == collapsed_bkd_path_lbndry[m, 2]]$uncollapsed
        cnt_bck_hap_b <- length(idx_set_hap_b)
        smp_prob <- as.vector(marginal_uncollapsed_fwd_prob_lbndry[idx_set_hap_b] / marginal_collapsed_fwd_prob_lbndry[collapsed_bkd_path_lbndry[m, 2]])
        uncollapsed_bkd_path_lbndry[m, 2] <- idx_set_hap_b[sample.int(cnt_bck_hap_b, size = 1, replace = TRUE, prob = smp_prob)]
      }
      
      return(uncollapsed_bkd_path_lbndry)
    }
    cmpuncollapseBkdPath_Diplo_FullRed_Bck <- cmpfun(uncollapseBkdPath_Diplo_FullRed_Bck)
    
    collapsed_ref_bck <- collapsed_info$unq
    n_bck <- collapsed_info$cnt
    N_bck <- nrow(collapsed_ref_bck)
    L_bck <- ncol(collapsed_ref_bck)
    
    collapsed_bkd_path <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_bck))
    
    # print(paste("marker:", L_bck))
    collapsed_bkd_path[, , L_bck] <- cmpcollapseBkdPath_Diplo_FullRed_Bck(uncollapsed_bkd_path_rbndry, smp_dip_size, collapsed_info)
    # smp_dip[, , L_bck] <- matrix(NA, nrow = smp_dip_size, ncol = 2)
    
    for (l in L_bck:2) {
      # print(paste("marker:", l - 1))
      for (m in 1:smp_dip_size) {
        trans_prob <- cmpcalculateTransProb_Diplo_Red_Pnl(collapsed_bkd_path[m, 1, l], collapsed_bkd_path[m, 2, l], theta_bck[l - 1], n_bck, N_bck, N_pnl)
        smp_prob <- as.vector(collapsed_fwd_prob[, , l - 1] * trans_prob / sum(collapsed_fwd_prob[, , l - 1] * trans_prob))
        collapsed_bkd_path_idx <- sample.int(N_bck * N_bck, size = 1, replace = TRUE, prob = smp_prob)
        collapsed_bkd_path[m, , l - 1] <- 
          if (collapsed_bkd_path_idx %% N_bck == 0) c(N_bck, collapsed_bkd_path_idx %/% N_bck) else c(collapsed_bkd_path_idx %% N_bck, collapsed_bkd_path_idx %/% N_bck + 1)
        smp_dip[m, , l - 1] <- cmpgenerateDiplotype_Diplo(obs_bck[l - 1], collapsed_ref_bck[collapsed_bkd_path[m, , l - 1], l - 1], lambda)
      }
      # print(smp_dip[, , l - 1])
    }
    uncollapsed_bkd_path_lbndry <- cmpuncollapseBkdPath_Diplo_FullRed_Bck(collapsed_bkd_path[, , 1], 
                                                                          collapsed_fwd_prob[, , 1], marginal_uncollapsed_fwd_prob_lbndry, smp_dip_size, collapsed_info)
    
    return(list(diplotype = smp_dip[, , 1:(L_bck - 1)], uncollapsed_bkd_path_lbndry = uncollapsed_bkd_path_lbndry))
  }
  cmprunBkdProcedure_Diplo_Red_Bck <- cmpfun(runBkdProcedure_Diplo_Red_Bck)
  
  # run the backward procedure
  uncollapsed_bkd_path_bndry <- array(NA, dim = c(smp_dip_size, 2, K_bck + 1))
  smp_dip <- array(NA, dim = c(smp_dip_size, 2, L_pnl))
  
  smp_prob <- as.vector(uncollapsed_fwd_prob_bndry[, K_bck + 1] / sum(uncollapsed_fwd_prob_bndry[, K_bck + 1]))
  uncollapsed_bkd_path_bndry[, 1, K_bck + 1] <- sample.int(N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  uncollapsed_bkd_path_bndry[, 2, K_bck + 1] <- sample.int(N_pnl, size = smp_dip_size, replace = TRUE, prob = smp_prob)
  for (m in 1:smp_dip_size) {
    smp_dip[m, , L_pnl] <- 
      cmpgenerateDiplotype_Diplo(obs_gen_ptn[[K_bck]][ncol(ref_pnl_ptn[[K_bck]])], ref_pnl_ptn[[K_bck]][uncollapsed_bkd_path_bndry[m, , K_bck + 1], ncol(ref_pnl_ptn[[K_bck]])], lambda)
  }
  
  bck_bndry <- L_pnl
  for (k in K_bck:1) {
    print(paste("block:", k))
    L_bck <- ncol(ref_pnl_ptn[[k]])
    collapsed_info <- cmpcollapseRefBlock(ref_pnl_ptn[[k]])
    bkd_procedure_bck <- cmprunBkdProcedure_Diplo_Red_Bck(obs_gen_ptn[[k]], ref_pnl_ptn[[k]], collapsed_info, theta_ptn[[k]], lambda, N_pnl, 
                                                          uncollapsed_bkd_path_bndry[, , k + 1], collapsed_fwd_prob[[k]], uncollapsed_fwd_prob_bndry[, k], smp_dip_size)
    smp_dip[, , (bck_bndry - L_bck  + 1):(bck_bndry - 1)] <- bkd_procedure_bck$diplotype
    uncollapsed_bkd_path_bndry[, , k] <- bkd_procedure_bck$uncollapsed_bkd_path_lbndry
    bck_bndry <- bck_bndry - L_bck  + 1
  }
  
  return(smp_dip)
}
#' Compiled version
cmprunBkdProcedure_Diplo_Red_Pnl_TA <- cmpfun(runBkdProcedure_Diplo_Red_Pnl_TA)

#########################

#' Run the block-wise path sampling with state-space reduction using Eqs.13-29, 60-74, 85-86 and 89-92 (the third type of the approximation)
#' @param obs_smp the observed sample
#' @param ref_pnl the reference panel
#' @param bck_len the length of the genomic block
#' @param theta the switch rates across the chromosome
#' @param lambda the mutation rate
#' @param smp_dip_size the number of the sampled diplotypes
#' @return the estimated diplotypes across the chromosome achieving the maximum a posterior returned in an array (dim = c(sample, haplotype, marker))

#' Standard version
runBWPS_Diplo_Red_Pnl_TA <- function(obs_smp, ref_pnl, bck_len, theta, lambda, smp_dip_size) {
  message("BWPS with state-space reduction using the 3rd type of the approximation at the boundary of each block")
  N_smp <- nrow(obs_smp)
  L_smp <- ncol(obs_smp)
  
  # run the preprocessing
  ref_pnl_ptn <- cmppartitionRefPanel(ref_pnl, bck_len, bck_idx = NULL)
  theta_ptn <- cmppartitionSwitchRate(theta, bck_len, bck_idx = NULL)
  
  N_pnl <- nrow(ref_pnl)
  L_pnl <- ncol(ref_pnl)
  K_bck <- length(ref_pnl_ptn)
  
  dip <- array(NA, dim = c(N_smp, 2, L_smp))
  for (i in 1:N_smp) {
    print(paste("individual:", i))
    obs_gen_ptn <- cmppartitionObsGenotype(obs_smp[i, ], bck_len, bck_idx = NULL)
    fwd_prob <- cmprunFwdProcedure_Diplo_Red_Pnl_TA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, K_bck)
    smp_dip <- cmprunBkdProcedure_Diplo_Red_Pnl_TA(obs_gen_ptn, ref_pnl_ptn, theta_ptn, lambda, N_pnl, L_pnl, K_bck, 
                                                   fwd_prob$collapsed_fwd_prob, fwd_prob$uncollapsed_fwd_prob_bndry, smp_dip_size)
    smp_dip_dist <- cmpconstructHist_Diplo(smp_dip, smp_dip_size)
    dip[i, , ] <- smp_dip_dist$diplotype[which(smp_dip_dist$count == max(smp_dip_dist$count))[1], , ]
    
    save(dip, file = "BWPS_Diplo_Red_Pnl_TA_N1000_L0500.rda")
  }
  return(dip)
}
#' Compiled version
cmprunBWPS_Diplo_Red_Pnl_TA <- cmpfun(runBWPS_Diplo_Red_Pnl_TA)

############################################################################################################################

#' @section Section 3 Results

#' Get the pattern of the forward probability matrix
#' @param fwd_prob the forward probability matrix
#' @return the pattern of the forward probability matrix

#' Standard version
getFwdProbPattern <- function(fwd_prob) {
  fwd_prob <- signif(fwd_prob, digits = 6)
  
  fwd_prob_pat <- fwd_prob
  fwd_prob_val <- sort(unique(as.vector(fwd_prob)))
  for (i in 1:length(fwd_prob_val)) {
    hap_a_idx <- which(fwd_prob == fwd_prob_val[i], arr.ind = TRUE)[, 1]
    hap_b_idx <- which(fwd_prob == fwd_prob_val[i], arr.ind = TRUE)[, 2]
    for (j in 1:length(hap_a_idx)) {
      fwd_prob_pat[hap_a_idx[j], hap_b_idx[j]] <- i
    }
  }
  
  return(fwd_prob_pat)
}
#' Compiled version
cmpgetFwdProbPattern <- cmpfun(getFwdProbPattern)

#########################

#' Calculate the switch error
#' @param true_dip the true diplotype
#' @param est_dip the estimated diplotype
#' @return the switch error including the markers required to switch returned in a list

#' Standard version
calculateSwitchError <- function(true_dip, est_dip) {
  calculateSwitchError_Haplo <- function(true_hap, est_hap) {
    hap_len <- length(true_hap)
    
    difference_SNP <- true_hap == est_hap
    
    switch_err <- 0
    switch_mkr <- rep(0, hap_len)
    
    if (hap_len == 1) {
      if (difference_SNP[1] == TRUE) {
        switch_err <- switch_err
        switch_mkr[1] <- 0
      } else {
        switch_err <- switch_err + 1
        switch_mkr[1] <- 1
      }
    } else {
      if (difference_SNP[1] == TRUE) {
        switch_err <- switch_err
        switch_mkr[1] <- 0
      } else {
        switch_err <- switch_err + 1
        switch_mkr[1] <- 1
      }
      
      for (l in 2:hap_len) {
        if (difference_SNP[l] == difference_SNP[l - 1]) {
          switch_err <- switch_err
          switch_mkr[l] <- 0
        } else {
          switch_err <- switch_err + 1
          switch_mkr[l] <- 1
        }
      }
    }
    
    return(list(switch_error = switch_err, switch_marker = switch_mkr))
  }
  cmpcalculateSwitchError_Haplo <- cmpfun(calculateSwitchError_Haplo)
  
  hetero_mkr <- which(colSums(true_dip) == 1)
  switch_err_case_1 <- cmpcalculateSwitchError_Haplo(true_dip[1, hetero_mkr], est_dip[1, hetero_mkr])
  switch_err_case_2 <- cmpcalculateSwitchError_Haplo(true_dip[1, hetero_mkr], est_dip[2, hetero_mkr])
  if (switch_err_case_1$switch_error < switch_err_case_2$switch_error) {
    switch_mkr <- rep(0, ncol(true_dip))
    switch_mkr[hetero_mkr] <- switch_err_case_1$switch_marker
    
    return(list(switch_error = switch_err_case_1$switch_error, switch_marker = switch_mkr))
  } else {
    switch_mkr <- rep(0, ncol(true_dip))
    switch_mkr[hetero_mkr] <- switch_err_case_2$switch_marker
    
    return(list(switch_error = switch_err_case_2$switch_error, switch_marker = switch_mkr))
  }
}
#' Compiled version
cmpcalculateSwitchError <- cmpfun(calculateSwitchError)

############################################################################################################################

#' @section Section 4 Discussion

############################################################################################################################

#' @section Appendix

