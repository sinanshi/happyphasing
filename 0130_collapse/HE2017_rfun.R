#' @title An accurate and efficient statistical method for phasing haplotypes using large phased reference panels
#' @author Zhangyi He

#' R functions

#install.packages("data.table")
library("data.table")

#install.packages("inline")
#library("inline")
#install.packages("Rcpp")
#library("Rcpp")

#install.packages("compiler")
#library("compiler")
#enableJIT(1)

#' call C++ functions
#sourceCpp("./Code/Code v1/HE2017_cfun.cpp")

#############################################################################################################################

#' @section Section 1. Introduction

#############################################################################################################################

#' @section Section 2. Materials and Methods

#' @section Section 2.1. Hidden Markov model

#' @section Section 2.2. Haplotype phasing with state-space reduction

####################################################################################################

#' Initialisation
#' Generate a haplotype reference panel
#' @param N_panel the number of the haplotypes in the reference panel
#' @param L_panel the number of the sites across the haplotype in the reference panel
#' @return a haplotype reference panel returned in a matrix (row = haplotypes and col = sites)
generateRefPanel <- function(N_panel, L_panel, seed = 1) {
  set.seed(seed)
  
  ref_panel <- matrix(as.integer(sample(c(0, 1), N_panel * L_panel, TRUE, rep(0.5, 2))), N_panel, L_panel)  
  
  return(ref_panel)
}

#' Generate a observed genotype sample
#' @param N_sample the number of the individuals (genotypes) in the observed sample
#' @param L_sample the number of the sites across the genotype in the observed sample
#' @return a observed genotype sample returned in a matrix (row = genotypes and col = sites)
#' genotype ranges from 0 to 3 and 0 = rr, 1 = ra/ar, 2 = aa, 3 = ?
generateObsSample <- function(N_sample, L_sample, seed = 2) {
  set.seed(seed)
  
  obs_sample <- matrix(as.integer(sample(c(0, 1, 2, 3), N_sample * L_sample, TRUE, rep(0.25, 4))), N_sample, L_sample)
  
  return(obs_sample)
}

#' Generate the switch rates across the genome
#' @param N_panel the number of the haplotypes in the reference panel
#' @param Ne the effective population size
#' @param d the genetic distances between successive sites
#' @return the switch rates across the genome returned in a vector
generateSwitchRate <- function(N_panel, Ne, d) {
  rho <- 4 * Ne * d
  theta <- 1 - exp(-rho / N_panel)
  
  return(theta)
}

#' Generate the switch rate across the chromosome 20 with the HRC
#' @param ref_panel the reference panel for chromosome 20 created by the HRC
#' @param N_panel the number of the haplotypes in the reference panel
#' @param Ne the effective population size
#' @return the switch rate across the chromosome 20 with the HRS returned in a vector
generateSwitchRate_HRC_Chr20 <- function(ref_panel, N_panel, Ne) {
  genetic_map <- read.table("genetic_map_chr20_combined_b37.txt", header = TRUE)
  genetic_map <- approx(as.vector(genetic_map[, 1]), as.vector(genetic_map[, 3]), xout = as.numeric(as.vector(colnames(ref_panel))))
  genetic_map_cm <- genetic_map$y
  genetic_map_cm[which(genetic_map_cm < 0)] <- 0
  d <- diff(genetic_map_cm, lag = 1) / 100
  rho <- 4 * Ne * d
  theta <- 1 - exp(-rho / N_panel)
  
  return(theta)
}

#' Generate the mutation rate
#' @param N_panel the number of the haplotypes in the reference panel
#' @return the mutation rate
generateLambda <- function(N_panel) {
  theta <- 1 / sum(1 / (1:(N_panel - 1)))
  lambda <- 0.5 * theta / (N_panel + theta)
  
  return(lambda)
}

#' Generate the emission probability table with Tab.1 
#' @param lambda the mutation rate
#' @return the emission probability table returned in a matrix (row = hidden genotypes, col = observed genotypes)
generateEmnProbTab <- function(lambda) {
  emn_prob <- matrix(NA, 3, 4)
  
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

###########################################################################

#' Preprocessing

#' Partition the haplotype reference panel into several genomic blocks with a predetermined block length
#' @param ref_panel the reference panel
#' @param block_length the length of the genomic block
#' @param block_index the indices of the genomic block
#' @return the haplotype reference blocks returned in a list
partitionRefPanel <- function(ref_panel, block_length, block_index = NULL) {
  L_panel <- ncol(ref_panel)
  
  if (is.null(block_index)) {
    block_index <- seq(1, L_panel, block_length - 1)
    if (block_index[length(block_index)] < L_panel) {
      block_index <- c(block_index, L_panel)
    }
  }
  
  ref_block <- list()
  for (i in seq(length(block_index) - 1)) {
    ref_block[[length(ref_block) + 1]] <- ref_panel[, block_index[i]:block_index[i + 1]]
  }
  
  return(lapply(ref_block, as.matrix))
}

#' Partition the observed genotype individual into several genomic blocks with a predetermined block length
#' @param obs_individual the observed individual
#' @param block_length the length of the genomic block
#' @param block_index the indices of the genomic block
#' @return the observed genotype blocks returned in a list
partitionObsIndividual <- function(obs_individual, block_length, block_index = NULL) {
  L_individual <- length(obs_individual)
  
  if (is.null(block_index)) {
    block_index <- seq(1, L_individual, block_length - 1)
    if (block_index[length(block_index)] < L_individual) {
      block_index <- c(block_index, L_individual)
    }
  }
  
  obs_block <- list()
  for (i in seq(length(block_index) - 1)) {
    obs_block[[length(obs_block) + 1]] <- obs_individual[block_index[i]:block_index[i + 1]]
  }
  
  return(lapply(obs_block, as.vector))
}

#' Partition the switch rates across the genome into several genomic blocks with a predetermined block length
#' @param theta the switch rates across the genome
#' @param block_length the length of the genomic block
#' @param block_index the indices of the genomic block
#' @return the switch rate blocks returned in a list
partitionSwitchRate <- function(theta, block_length, block_index = NULL) {
  L_panel <- length(theta) + 1
  
  if (is.null(block_index)) {
    block_index <- seq(1, L_panel, block_length - 1)
    if (block_index[length(block_index)] < L_panel) {
      block_index <- c(block_index, L_panel)
    }
  }
  
  theta_block <- list()
  for (i in seq(length(block_index) - 1)) {
    theta_block[[length(theta_block) + 1]] <- theta[block_index[i]:(block_index[i + 1] - 1)]
  }
  
  return(lapply(theta_block, as.vector))
}

###########################################################################

#' @section Section 2.2.1 Forward procedure

#' Forward procedure for a diplotype case

#' Calculate the emission probabilities at a specific site with Tab.1
#' @param obs_genotype the genotype at a specific site of the observed individual
#' @param ref_haplotype1 the first haplotype at a specific site of the reference panel
#' @param ref_haplotype2 the second haplotype at a specific site of the reference panel
#' @param emn_prob_tab the emission probability table
#' @return the emission probabilities at a specific site returned as a matrix
calculateEmnProb_Diplo <- function(obs_genotype, ref_haplotype1, ref_haplotype2, emn_prob_tab) {
  n1 <- length(ref_haplotype1)
  n2 <- length(ref_haplotype2)
  
  ref_h1 <- matrix(ref_haplotype1, n1, n2)
  ref_h2 <- t(matrix(ref_haplotype2, n2, n1))
  ref_g <- ref_h1 + ref_h2
  
  emn_prob <- function(ref_genotype)
    return(emn_prob_tab[(ref_genotype + 1), (obs_genotype + 1)])
  
  return(apply(ref_g, 2, emn_prob))
}

##################################################

#' Forward procedure without state-space reduction

#' Perform the forward procedure using the haplotype reference panel without state-space reduction with Eqs.9 and 10
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param emn_prob_tab the emission probability table
#' @return the forward probabilities across the whole genome returned as a list
performFwdProcedure_Diplo_Std_Panel <- function(obs_individual, ref_panel, theta, emn_prob_tab) {
  N_panel <- nrow(ref_panel)
  L_panel <- ncol(ref_panel)
  
  fwd_prob <- list()
  
  # initialise the forward probabilities at the first site of the haplotype reference panel with Eq.10
  fwd_prob[[1]] <- calculateEmnProb_Diplo(obs_individual[1], ref_panel[, 1], ref_panel[, 1], emn_prob_tab) / (N_panel ^ 2)
  
  # update the forward probabilities across the haplotype reference panel with Eq.9
  for (l in 2:L_panel) {
    fwd_prob[[l]] <- matrix(NA, N_panel, N_panel)
    emn_prob <- calculateEmnProb_Diplo(obs_individual[l], ref_panel[, l], ref_panel[, l], emn_prob_tab)
    for (i in seq(N_panel)) {
      for (j in seq(N_panel)) {
        fwd_prob[[l]][i, j] <- 
          emn_prob[i, j] * ((1 - theta[l - 1]) ^ 2) * fwd_prob[[l - 1]][i, j] +
          emn_prob[i, j] * (1 - theta[l - 1]) * (theta[l - 1] / N_panel) * sum(fwd_prob[[l - 1]][i, ]) +
          emn_prob[i, j] * (theta[l - 1] / N_panel) * (1 - theta[l - 1]) * sum(fwd_prob[[l - 1]][, j]) +
          emn_prob[i, j] * (theta[l - 1] / N_panel) ^ 2 * sum(fwd_prob[[l - 1]])
      }
    }
  }
  
  return(fwd_prob)
}

#########################

#' Perform the forward procedure using the haplotype reference block without state-space reduction (including the forward probabilities in each case)

#' Initialise the forward probabilities at the second site of the genomic block
#' @param fwd_prob_lbndry the forward probabilities at the first site of the genomic block
#' @param emn_prob the emission probabilities at the second site of the genomic block
#' @param theta the switch rate between the first and second sites of the genomic block
#' @param N_block the number of the haplotypes in the reference block
#' @return the forward probabilities at the second site of the genomic block returned as a list
initialiseFwdProb_Diplo_Std_Block <- function(fwd_prob_lbndry, emn_prob, theta, N_block) {
  # case 1: non-recombination & non-recombination
  fwd_prob_nr_nr <- function() 
    emn_prob * (1 - theta) ^ 2 * fwd_prob_lbndry
  
  # case 2: non-recombination & recombination
  fwd_prob_nr_r <- function() 
    apply(emn_prob * (1 - theta) * (theta / N_block), 2, 
          function(x) x * rowSums(fwd_prob_lbndry))
  
  # case 3: recombination & non-recombination
  fwd_prob_r_nr <- function()
    apply(emn_prob %*% diag((theta / N_block) * (1 - theta) * colSums(fwd_prob_lbndry)), 2, 
          function(x) x)
  
  # case 4: recombination & recombination
  fwd_prob_r_r <- function()
    apply(emn_prob * (theta / N_block) * (theta / N_block) * sum(fwd_prob_lbndry), 2, 
          function(x) x)
  
  fwd_prob <- list()
  fwd_prob[["nr_nr"]] <- fwd_prob_nr_nr()
  fwd_prob[["nr_r"]] <- fwd_prob_nr_r()
  fwd_prob[["r_nr"]] <- fwd_prob_r_nr()
  fwd_prob[["r_r"]] <- fwd_prob_r_r()
  fwd_prob[["total"]] <- Reduce('+', fwd_prob)
  
  return(fwd_prob)
}

#' Update the forward probabilities at the current site of the genomic block
#' @param fwd_prob_pre the forward probabilities at the preceding site of the genomic block
#' @param emn_prob the emission probabilities at the current site of the genomic block
#' @param theta the switch rate between the preceding and current sites of the genomic block
#' @param N_block the number of the haplotypes in the reference block
#' @return the forward probabilities at the current site of the genomic block returned as a list
updateFwdProb_Diplo_Std_Block <- function(fwd_prob_pre, emn_prob, theta, N_block) {
  # case 1: non-recombination & non-recombination
  fwd_prob_nr_nr <- function() {
    emn_prob * (1 - theta) ^ 2 * fwd_prob_pre$nr_nr
  }
  
  # case 2: non-recombination & recombination
  fwd_prob_nr_r <- function() {
    emn_prob * (1 - theta) ^ 2 * fwd_prob_pre$nr_r + 
      apply(emn_prob * (1 - theta) * (theta / N_block), 2, 
            function(x) x * rowSums(fwd_prob_pre$nr_nr + fwd_prob_pre$nr_r))
  }
  
  # case 3: recombination & non-recombination
  fwd_prob_r_nr <- function() {
    emn_prob * (1 - theta) ^ 2 * fwd_prob_pre$r_nr + 
      apply(emn_prob %*% diag((theta / N_block) * (1 - theta) * colSums(fwd_prob_pre$nr_nr + fwd_prob_pre$r_nr)), 2, 
            function(x) x)
  }
  
  # case 4: recombination & recombination
  fwd_prob_r_r <- function() {
    emn_prob * (1 - theta) ^ 2 * fwd_prob_pre$r_r + 
      apply(emn_prob * (1 - theta) * (theta / N_block), 2, 
            function(x) x * rowSums(fwd_prob_pre$r_nr + fwd_prob_pre$r_r)) + 
      apply(emn_prob %*% diag((theta / N_block) * (1 - theta) * colSums(fwd_prob_pre$nr_r + fwd_prob_pre$r_r)), 2, 
            function(x) x) + 
      apply(emn_prob * (theta / N_block) * (theta / N_block) * sum(fwd_prob_pre$total), 2, 
            function(x) x)
  }
  
  fwd_prob <- list()
  fwd_prob[["nr_nr"]] <- fwd_prob_nr_nr()
  fwd_prob[["nr_r"]] <- fwd_prob_nr_r()
  fwd_prob[["r_nr"]] <- fwd_prob_r_nr()
  fwd_prob[["r_r"]] <- fwd_prob_r_r()
  fwd_prob[["total"]] <- Reduce('+', fwd_prob)
  
  return(fwd_prob)
}

#' Perform the forward procedure using the haplotype reference block without state-space reduction
#' @param obs_block the observed block
#' @param ref_block the reference block
#' @param theta the switch rates across the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the forward probabilities across the genomic block returned as a list
performFwdProcedure_Diplo_Std_Block <- function(fwd_prob_lbndry, obs_block, ref_block, theta, emn_prob_tab) {
  N_block <- nrow(ref_block)
  L_block <- ncol(ref_block)
  
  fwd_prob <- list()
  
  # initialise the forward probabilities at the first site of the genomic block
  fwd_prob[[1]] <- list(total = fwd_prob_lbndry)
  
  # initialise the forward probabilities at the second site of the genomic block
  emn_prob <- calculateEmnProb_Diplo(obs_block[2], ref_block[, 2], ref_block[, 2], emn_prob_tab)
  fwd_prob[[2]] <- initialiseFwdProb_Diplo_Std_Block(fwd_prob[[1]]$total, emn_prob, theta[1], N_block)
  
  # update the forward probabilities across the genomic block (site >= 3)
  for (l in 3:L_block) {
    emn_prob <- calculateEmnProb_Diplo(obs_block[l], ref_block[, l], ref_block[, l], emn_prob_tab)
    fwd_prob[[l]] <- updateFwdProb_Diplo_Std_Block(fwd_prob[[l - 1]], emn_prob, theta[l - 1], N_block)
  }
  
  return(fwd_prob)
}

##################################################

#' Forward procedure with state-space reduction

#' Collapse the haplotype reference block (i.e., remove the duplicates of the haplotypes in the haplotype reference block)
#' @param ref_block the reference block
#' @return the collapsed haplotype reference block with the associated map to the uncollapsed haplotype reference block returned in a list
collapseRefBlock <- function(ref_block) {
  unq <- unique(ref_block)
  
  env <- new.env()
  env$ind <- list()
  
  cnt <- vector()
  cnt <- apply(unq, 1, 
               function(x) {
                 ind_collapsed <- length(env$ind) + 1
                 ind_uncollapsed <- which(apply((t(t(ref_block) - x)) == 0, 1, all))
                 env$ind[[ind_collapsed]] <- ind_uncollapsed
                 env$map[[ind_collapsed]] <- data.frame(collapsed = rep(ind_collapsed, length(ind_uncollapsed)), uncollapsed = ind_uncollapsed)
                 return(length(ind_uncollapsed))
               }
  )
  
  map <- rbindlist(env$map)
  setkey(map, "uncollapsed")
  
  return(list(unq = unq, cnt = cnt, ind = env$ind, map = map))
}

#' Initialise the forward probabilities at the first site of the whole genome with Eq.10
#' @param obs_genotype the genotype at the first site of the observed individual
#' @param ref_haplotype the haplotypes at the first site of the reference panel
#' @param N_panel the number of the haplotypes in the reference panel
#' @param emn_prob_tab the emission probability table
#' @return the forward probabilities at the first site of the whole genome returned as a matrix
initialiseFwdProb_Diplo_Red_Panel <- function(obs_genotype, ref_haplotype, N_panel, emn_prob_tab) {
  fwd_prob <- calculateEmnProb_Diplo(obs_genotype, ref_haplotype, ref_haplotype, emn_prob_tab) / (N_panel ^ 2)
  
  return(fwd_prob)
}

#########################

#' Perform the forward procedure using the haplotype reference block with the full-collapsed state space for the non-recombination & non-recombination case and the recombination & recombination case

#' Collapse the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block with Eq.16
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block returned as a matrix
collapseFwdProb_Diplo_FullRed_Block <- function(uncollapsed_fwd_prob_lbndry, collapsed_info) {
  N_block <- length(collapsed_info$cnt)
  
  collapsed_fwd_prob_lbndry <- matrix(NA, N_block, N_block)
  for (u in seq(N_block)) {
    for (v in seq(N_block)) {
      collapsed_fwd_prob_lbndry[u, v] <- sum(uncollapsed_fwd_prob_lbndry[collapsed_info$ind[[u]], collapsed_info$ind[[v]]])
    }
  }
  
  return(collapsed_fwd_prob_lbndry)
}

#' Initialise the full-collapsed forward probabilities at the second site of the genomic block with Eqs.14, 22, 24, 26 and 28
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob the emission probabilities at the second site of the genomic block
#' @param theta the switch rate between the first and second sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @return the full-collapsed forward probabilities at the second site of the genomic block returned as a list
initialiseFwdProb_Diplo_FullRed_Block <- function(collapsed_fwd_prob_lbndry, emn_prob, theta, N_panel, N_block) {
  # case 1: non-recombination & non-recombination
  collapsed_fwd_prob_nr_nr <- function() 
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_lbndry
  
  # case 2: non-recombination & recombination
  collapsed_fwd_prob_nr_r <- function() 
    apply(emn_prob %*% diag((1 - theta) * (theta * N_block / N_panel)), 2, 
          function(x) x * rowSums(collapsed_fwd_prob_lbndry))
  
  # case 3: recombination & non-recombination
  collapsed_fwd_prob_r_nr <- function()
    apply(emn_prob %*% diag((theta / N_panel) * (1 - theta) * colSums(collapsed_fwd_prob_lbndry)), 2, 
          function(x) x * N_block)
  
  # case 4: recombination & recombination
  collapsed_fwd_prob_r_r <- function()
    apply(emn_prob %*% diag((theta * N_block / N_panel) * (theta / N_panel) * sum(collapsed_fwd_prob_lbndry)), 2, 
          function(x) x * N_block)
  
  collapsed_fwd_prob <- list()
  collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
  collapsed_fwd_prob[["nr_r"]] <- collapsed_fwd_prob_nr_r()
  collapsed_fwd_prob[["r_nr"]] <- collapsed_fwd_prob_r_nr()
  collapsed_fwd_prob[["r_r"]] <- collapsed_fwd_prob_r_r()
  collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
  
  return(collapsed_fwd_prob)
}

#' Update the full-collapsed forward probabilities at the current site of the genomic block with Eqs.14, 21, 23, 25 and 27
#' @param collapsed_fwd_prob_pre the full-collapsed forward probabilities at the preceding site of the genomic block
#' @param emn_prob the emission probabilities at the current site of the genomic block 
#' @param theta the switch rate between the preceding and current sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @return the full-collapsed forward probabilities at the current site of the genomic block returned as a list
updateFwdProb_Diplo_FullRed_Block <- function(collapsed_fwd_prob_pre, emn_prob, theta, N_panel, N_block) {
  # case 1: non-recombination & non-recombination
  collapsed_fwd_prob_nr_nr <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
  }
  
  # case 2: non-recombination & recombination  
  collapsed_fwd_prob_nr_r <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_r + 
      apply(emn_prob %*% diag((1 - theta) * (theta * N_block / N_panel)), 2, 
            function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_r))
  }
  
  # case 3: recombination & non-recombination
  collapsed_fwd_prob_r_nr <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$r_nr + 
      apply(emn_prob %*% diag((theta / N_panel) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$r_nr)), 2, 
            function(x) x * N_block)
  }
  
  # case 4: recombination & recombination
  collapsed_fwd_prob_r_r <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$r_r + 
      apply(emn_prob %*% diag((1 - theta) * (theta * N_block / N_panel)), 2, 
            function(x) x * rowSums(collapsed_fwd_prob_pre$r_nr + collapsed_fwd_prob_pre$r_r)) + 
      apply(emn_prob %*% diag((theta / N_panel) * (1 - theta) * colSums(collapsed_fwd_prob_pre$nr_r + collapsed_fwd_prob_pre$r_r)), 2, 
            function(x) x * N_block) + 
      apply(emn_prob %*% diag((theta * N_block / N_panel) * (theta / N_panel) * sum(collapsed_fwd_prob_pre$total)), 2, 
            function(x) x * N_block)
  }
  
  collapsed_fwd_prob <- list()
  collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
  collapsed_fwd_prob[["nr_r"]] <- collapsed_fwd_prob_nr_r()
  collapsed_fwd_prob[["r_nr"]] <- collapsed_fwd_prob_r_nr()
  collapsed_fwd_prob[["r_r"]] <- collapsed_fwd_prob_r_r()
  collapsed_fwd_prob[["total"]] <- Reduce('+', collapsed_fwd_prob)
  
  return(collapsed_fwd_prob)
}

#' Uncollapse the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block for the non-recombination & non-recombination case and the recombination & recombination case with Eqs.30 and 31
#' @param collapsed_fwd_prob_rbndry the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & non-recombination case and the recombination & recombination case returned as a list
uncollapseFwdProb_Diplo_FullRed_Block <- function(collapsed_fwd_prob_rbndry, collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info) {
  uncollapsed_fwd_prob_rbndry_nr_nr <- matrix(NA, N_panel, N_panel)
  uncollapsed_fwd_prob_rbndry_r_r <- matrix(NA, N_panel, N_panel)
  for (i in 1:N_panel) {
    for (j in 1:N_panel) {
      u <- collapsed_info$map[uncollapsed == i]$collapsed
      v <- collapsed_info$map[uncollapsed == j]$collapsed
      N_block_u <- collapsed_info$cnt[u]
      N_block_v <- collapsed_info$cnt[v]
      uncollapsed_fwd_prob_rbndry_nr_nr[i, j] <- collapsed_fwd_prob_rbndry$nr_nr[u, v] * uncollapsed_fwd_prob_lbndry[i, j] / collapsed_fwd_prob_lbndry[u, v]
      uncollapsed_fwd_prob_rbndry_r_r[i, j] <- collapsed_fwd_prob_rbndry$r_r[u, v] / N_block_u / N_block_v
    }
  }
  
  return(list(nr_nr = uncollapsed_fwd_prob_rbndry_nr_nr, r_r = uncollapsed_fwd_prob_rbndry_r_r))
}

#' Perform the forward procedure using the haplotype reference block with the full-collapsed state space
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the full-collapsed forward probabilities across the genomic block and the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block for the non-recombination & non-recombination case and the recombination & recombination case returned as a list
performFwdProcedure_Diplo_FullRed_Block <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_fwd_prob <- list()
  
  # calculate the full-collapsed forward probabilities at the first site of the genomic block
  collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
  
  if (L_block > 1) {
    # calculate the full-collapsed forward probabilities at the second site of the genomic block
    collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
    collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
    
    if (L_block > 2) {
      # calculate the full-collapsed forward probabilities at the following sites of the genomic block
      for (l in 3:L_block) {
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
        collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
      }
    }
  }
  
  #' calculate the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & non-recombination case and the recombination & recombination case
  uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
  uncollapsed_fwd_prob_rbndry_nr_nr <- uncollapsed_fwd_prob_rbndry$nr_nr
  uncollapsed_fwd_prob_rbndry_r_r <- uncollapsed_fwd_prob_rbndry$r_r
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry_nr_nr = uncollapsed_fwd_prob_rbndry_nr_nr, uncollapsed_fwd_prob_rbndry_r_r = uncollapsed_fwd_prob_rbndry_r_r))
}

#########################

#' Perform the forward procedure using the haplotype reference block with the semi-collapsed state space for the non-recombination & recombination case (or the recombination & non-recombination case)

#' Collapse the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block with Eq.35 (or 48)
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the semi-collapsed forward probabilities at the first site (the left boundary) of the genomic block returned as a matrix
collapseFwdProb_Diplo_SemiRed_Block <- function(uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info) {
  N_block <- length(collapsed_info$cnt)
  
  collapsed_fwd_prob_lbndry <- matrix(NA, N_panel, N_block)
  for (i in seq(N_panel)) {
    for (v in seq(N_block)) {
      collapsed_fwd_prob_lbndry[i, v] <- sum(uncollapsed_fwd_prob_lbndry[i, collapsed_info$ind[[v]]])
    }
  }
  
  return(collapsed_fwd_prob_lbndry)
}

#' Initialise the semi-collapsed forward probabilities at the second site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case with Eqs.41 and 43 (or 54 and 56)
#' @param collapsed_fwd_prob_lbndry the semi-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob the emission probabilities at the second site of the genomic block
#' @param theta the switch rate between the first and second sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @return the semi-collapsed forward probabilities at the second site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case returned as a list
initialiseFwdProb_Diplo_SemiRed_Block <- function(collapsed_fwd_prob_lbndry, emn_prob, theta, N_panel, N_block) {
  # case 1: non-recombination & non-recombination
  collapsed_fwd_prob_nr_nr <- function() 
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_lbndry
  
  # case 2: non-recombination & recombination
  collapsed_fwd_prob_nr_r <- function() 
    apply(emn_prob %*% diag((1 - theta) * (theta * N_block / N_panel)), 2, 
          function(x) x * rowSums(collapsed_fwd_prob_lbndry))
  
  collapsed_fwd_prob <- list()
  collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
  collapsed_fwd_prob[["nr_r"]] <- collapsed_fwd_prob_nr_r()
  
  return(collapsed_fwd_prob)
}

#' Update the semi-collapsed forward probabilities at the current site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case with Eqs.40 and 42 (53 and 55)
#' @param collapsed_fwd_prob_pre the semi-collapsed forward probabilities at the preceding site of the genomic block
#' @param emn_prob the emission probabilities at the current site of the genomic block 
#' @param theta the switch rate between the preceding and current sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @return the semi-collapsed forward probabilities at the current site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case returned as a list
updateFwdProb_Diplo_SemiRed_Block <- function(collapsed_fwd_prob_pre, emn_prob, theta, N_panel, N_block) {
  # case 1: non-recombination & non-recombination
  collapsed_fwd_prob_nr_nr <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_nr
  }
  
  # case 2: non-recombination & recombination
  collapsed_fwd_prob_nr_r <- function() {
    emn_prob * (1 - theta) ^ 2 * collapsed_fwd_prob_pre$nr_r + 
      apply(emn_prob %*% diag((1 - theta) * (theta * N_block / N_panel)), 2, 
            function(x) x * rowSums(collapsed_fwd_prob_pre$nr_nr + collapsed_fwd_prob_pre$nr_r))
  }
  
  collapsed_fwd_prob <- list()
  collapsed_fwd_prob[["nr_nr"]] <- collapsed_fwd_prob_nr_nr()
  collapsed_fwd_prob[["nr_r"]] <- collapsed_fwd_prob_nr_r()
  
  return(collapsed_fwd_prob)
}

#' Uncollapse the semi-collapsed forward probabilities at the last site (the right boundary) of the genomic block for the non-recombination & recombination case (or the recombination & non-recombination case) with Eq.44 (or 57)
#' @param collapsed_fwd_prob_rbndry the semi-collapsed forward probabilities at the last site (the right boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & recombination case (or the recombination & non-recombination case) returned as a list
uncollapseFwdProb_Diplo_SemiRed_Block <- function(collapsed_fwd_prob_rbndry, N_panel, collapsed_info) {
  uncollapsed_fwd_prob_rbndry_nr_r <- matrix(NA, N_panel, N_panel)
  for (i in 1:N_panel) {
    for (j in 1:N_panel) {
      v <- collapsed_info$map[uncollapsed == j]$collapsed
      N_block_v <- collapsed_info$cnt[v]
      uncollapsed_fwd_prob_rbndry_nr_r[i, j] <- collapsed_fwd_prob_rbndry$nr_r[i, v] / N_block_v
    }
  }
  
  return(list(nr_r = uncollapsed_fwd_prob_rbndry_nr_r, r_nr = t(uncollapsed_fwd_prob_rbndry_nr_r)))
}

#' Perform the forward procedure using the haplotype reference block with the semi-collapsed state space
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block for the non-recombination & recombination case or (the recombination & non-recombination case) returned as a list
performFwdProcedure_Diplo_SemiRed_Block <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_fwd_prob <- list()
  
  # calculate the semi-collapsed forward probabilities at the first site of the genomic block
  collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_SemiRed_Block(uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info))
  
  if (L_block > 1) {
    # calculate the semi-collapsed forward probabilities at the second site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
    collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], uncollapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
    collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_SemiRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
    
    if (L_block > 2) {
      # calculate the semi-collapsed forward probabilities at the following sites of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
      for (l in 3:L_block) {
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], uncollapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
        collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_SemiRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
      }
    }
  }
  
  #' calculate the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & recombination case
  uncollapsed_fwd_prob_rbndry_nr_r <- uncollapseFwdProb_Diplo_SemiRed_Block(collapsed_fwd_prob[[L_block]], N_panel, collapsed_info)$nr_r
  
  return(list(uncollapsed_fwd_prob_rbndry_nr_r = uncollapsed_fwd_prob_rbndry_nr_r))
}

#########################

#' Perform the forward procedure using the haplotype reference block with state-space reduction
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the genomic block and the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block returned as a list
performFwdProcedure_Diplo_Red_Block <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  # implement the HMM model over the collapsed state space
  full_collapsed_case <- performFwdProcedure_Diplo_FullRed_Block(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab)
  uncollapsed_fwd_prob_rbndry_nr_nr <- full_collapsed_case$uncollapsed_fwd_prob_rbndry_nr_nr
  uncollapsed_fwd_prob_rbndry_r_r <- full_collapsed_case$uncollapsed_fwd_prob_rbndry_r_r
  
  # implement the HMM model over the semi-collapsed state space
  semi_collapsed_case <- performFwdProcedure_Diplo_SemiRed_Block(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab)
  uncollapsed_fwd_prob_rbndry_nr_r <- semi_collapsed_case$uncollapsed_fwd_prob_rbndry_nr_r
  uncollapsed_fwd_prob_rbndry_r_nr <- t(semi_collapsed_case$uncollapsed_fwd_prob_rbndry_nr_r)
  
  collapsed_fwd_prob <- full_collapsed_case$collapsed_fwd_prob
  uncollapsed_fwd_prob_rbndry <- uncollapsed_fwd_prob_rbndry_nr_nr + uncollapsed_fwd_prob_rbndry_r_r + uncollapsed_fwd_prob_rbndry_nr_r + uncollapsed_fwd_prob_rbndry_r_nr
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
}

#' Perform the forward procedure using the haplotype reference panel with state-space reduction
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param block_length the length of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the whole genome and the uncollapsed forward probabilities at the boundaries of the genome blocks returned as a list
performFwdProcedure_Diplo_Red_Panel <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
  N_panel <- nrow(ref_panel)
  
  # execute the preprocessing
  obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
  ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
  theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
  
  K_panel <- length(ref_panel_partition)
  
  # execute the forward procedure
  uncollapsed_fwd_prob_bndry <- list()
  uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
  
  collapsed_fwd_prob <- list()
  for (k in 1:K_panel) {
    fwd_procedure_block <- performFwdProcedure_Diplo_Red_Block(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}

###########################################################################

#' @section Section 2.2.2 Backward procedure

#' Calculate the transition probabilities at a specific site with Eq.20
#' @param ref_haplo_index1 the first haplotype index at a specific site of the reference panel
#' @param ref_haplo_index2 the second haplotype index at a specific site of the reference panel
#' @param theta the switch rate between the preceding and current sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @return the transition probabilities at a specific site returned as a matrix
calculateTransProb_Diplo <- function(ref_haplo_index1, ref_haplo_index2, theta, N_panel, N_block) {
  trans_prob <- matrix((theta * N_block[ref_haplo_index1] / N_panel) * (theta * N_block[ref_haplo_index2] / N_panel), length(N_block), length(N_block))
  trans_prob[ref_haplo_index1, ] <- (theta * N_block[ref_haplo_index1] / N_panel + (1 - theta)) * (theta * N_block[ref_haplo_index2] / N_panel)
  trans_prob[, ref_haplo_index2] <- (theta * N_block[ref_haplo_index1] / N_panel) * (theta * N_block[ref_haplo_index2] / N_panel + (1 - theta))
  trans_prob[ref_haplo_index1, ref_haplo_index2] <- (theta * N_block[ref_haplo_index1] / N_panel + (1 - theta)) * (theta * N_block[ref_haplo_index2] / N_panel + (1 - theta))
  
  return(trans_prob)
}

#' Generate the diplotype at a specific site with Eqs.69-73
#' @param obs_genotype the genotype at a specific site of the observed individual
#' @param ref_haplotype1 the first haplotype at a specific site of the reference panel
#' @param ref_haplotype2 the second haplotype at a specific site of the reference panel
#' @param lambda the mutation rate
#' @return the diplotype at a specific site returned as a vector
generateDiplotype_Diplo <- function(obs_genotype, ref_haplotype1, ref_haplotype2, lambda) {
  diplotype <- numeric(2)
  
  if (obs_genotype == 0) {
    diplotype <- c(0, 0)
  }
  if (obs_genotype == 1) {
    if ((ref_haplotype1 + ref_haplotype2) == 1) {
      sampling_prob <- c((1 - lambda) ^ 2, lambda ^ 2) / ((1 - lambda) ^ 2 + lambda ^ 2)
      choice <- sample.int(2, size = 1, replace = TRUE, prob = sampling_prob)
      if (choice == 1) {
        diplotype <- c(ref_haplotype1, ref_haplotype2)
      }
      if (choice == 2) {
        diplotype <- c(1 - ref_haplotype1, 1 - ref_haplotype2)
      }
    }
    else {
      sampling_prob <- c(0.5, 0.5)
      choice <- sample.int(2, size = 1, replace = TRUE, prob = sampling_prob)
      if (choice == 1) {
        diplotype <- c(0, 1)
      }
      if (choice == 2) {
        diplotype <- c(1, 0)
      }
    }
  }
  if (obs_genotype == 2) {
    diplotype <- c(1, 1)
  }
  if (obs_genotype == 3) {
    sampling_prob <- c((1 - lambda) ^ 2, (1 - lambda) * lambda, lambda * (1 - lambda), lambda ^ 2)
    choice <- sample.int(4, size = 1, replace = TRUE, prob = sampling_prob)
    if (choice == 1) {
      diplotype <- c(ref_haplotype1, ref_haplotype2)
    }
    if (choice == 2) {
      diplotype <- c(ref_haplotype1, 1 - ref_haplotype2)
    }
    if (choice == 3) {
      diplotype <- c(1 - ref_haplotype1, ref_haplotype2)
    } 
    if (choice == 4) {
      diplotype <- c(1 - ref_haplotype1, 1 - ref_haplotype2)
    }  
  }
  
  return(diplotype)
}

##################################################

#' Backward procedure without state-space reduction

#' Perform the backward procedure using the haplotype reference panel without state-space reduction
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param fwd_prob the forward probabilities across the whole genome
#' @param theta the switch rates across the whole genome
#' @param lambda the mutation rate
#' @param M_sample the number of the backward paths generated
#' @return the backward path across the whole genome and the diplotype across the whole genome returned as a list
performBwdProcedure_Diplo_Std_Panel <- function(obs_individual, ref_panel, fwd_prob, theta, lambda, M_sample) {
  N_panel <- nrow(ref_panel)
  L_panel <- ncol(ref_panel)
  
  bwd_path <- list()
  diplotype <- array(NA, dim = c(2, L_panel, M_sample))
  
  # initialise the backward paths and the diplotypes at the last site of the haplotype reference panel
  bwd_path[[L_panel]] <- matrix(NA, 2, M_sample)
  for (m in 1:M_sample) {
    sampling_prob <- as.vector(fwd_prob[[L_panel]] / sum(fwd_prob[[L_panel]]))
    bwd_path_index <- sample.int(N_panel * N_panel, size = 1, replace = TRUE, prob = sampling_prob)
    bwd_path[[L_panel]][, m] <- if (bwd_path_index %% N_panel == 0) c(N_panel, bwd_path_index %/% N_panel) else c(bwd_path_index %% N_panel, bwd_path_index %/% N_panel + 1)
    ref_haplotype1 <- ref_panel[bwd_path[[L_panel]][1, m], L_panel]
    ref_haplotype2 <- ref_panel[bwd_path[[L_panel]][2, m], L_panel]
    diplotype[, L_panel, m] <- generateDiplotype_Diplo(obs_individual[L_panel], ref_haplotype1, ref_haplotype2, lambda)
  }
  
  # update the backward paths and the diplotypes across the haplotype reference panel
  if (L_panel > 1) {
    for (l in (L_panel - 1):1) {
      bwd_path[[l]] <- matrix(NA, 2, M_sample)
      for (m in 1:M_sample) {
        trans_prob <- calculateTransProb_Diplo(bwd_path[[l + 1]][1, m], bwd_path[[l + 1]][2, m], theta[l], N_panel, rep(1, N_panel))
        sampling_prob <- as.vector(fwd_prob[[l]] * trans_prob / sum(fwd_prob[[l]] * trans_prob))
        bwd_path_index <- sample.int(N_panel * N_panel, size = 1, replace = TRUE, prob = sampling_prob)
        bwd_path[[l]][, m] <- if (bwd_path_index %% N_panel == 0) c(N_panel, bwd_path_index %/% N_panel) else c(bwd_path_index %% N_panel, bwd_path_index %/% N_panel + 1)
        ref_haplotype1 <- ref_panel[bwd_path[[l]][1, m], l]
        ref_haplotype2 <- ref_panel[bwd_path[[l]][2, m], l]
        diplotype[, l, m] <- generateDiplotype_Diplo(obs_individual[l], ref_haplotype1, ref_haplotype2, lambda)
      }
    }
  }
  
  bwd_path <- aperm(array(unlist(bwd_path), dim = c(2, M_sample, L_panel)), c(1, 3, 2))
  
  return(list(bwd_path = bwd_path, diplotype = diplotype))
}

#########################

#' Perform the backward procedure using the haplotype reference block without state-space reduction

#' Update the backward path at the current site of the genomic block
#' @param
#' @return
#updateBwdProb_Diplo_Std_Block <- function()

#' Perform the backward path using the haplotype reference block without state-space reduction
#' @param
#' @return
#performBwdProcedure_Diplo_Std_Block <- function()

##################################################

#' Backward procedure with state-space reduction

#' Initialise the backward path at the first site of the whole genome with Eqs.60 and 61
#' @param obs_genotype the genotype at the first site of the observed individual
#' @param ref_haplotype the haplotypes at the last site of the reference panel
#' @param lambda the mutation rate
#' @param N_panel the number of the haplotypes in the reference panel
#' @param fwd_prob the forward probabilities at the last site of the whole genome
#' @param M_sample the number of the backward paths generated
#' @return the backward path and the diplotype at the last site of the whole genome returned as a list
initialiseBwdPath_Diplo_Red_Panel <- function(obs_genotype, ref_haplotype, lambda, N_panel, fwd_prob, M_sample) {
  bwd_path <- matrix(NA, 2, M_sample)
  diplotype <- matrix(NA, 2, M_sample)
  sampling_prob <- as.vector(fwd_prob / sum(fwd_prob))
  bwd_path_index <- sample.int(N_panel * N_panel, size = M_sample, replace = TRUE, prob = sampling_prob)
  for (m in 1:M_sample) {
    bwd_path[, m] <- if (bwd_path_index[m] %% N_panel == 0) c(N_panel, bwd_path_index[m] %/% N_panel) else c(bwd_path_index[m] %% N_panel, bwd_path_index[m] %/% N_panel + 1)
    ref_haplotype1 <- ref_haplotype[bwd_path[1, m]]
    ref_haplotype2 <- ref_haplotype[bwd_path[2, m]]
    diplotype[, m] <- generateDiplotype_Diplo(obs_genotype, ref_haplotype1, ref_haplotype2, lambda)
  }
  
  return(list(bwd_path = bwd_path, diplotype = diplotype))
}

#########################

#' Perform the backward procedure using the haplotype reference block with the full-collapsed state space

#' Collapse the uncollapsed backward path at the last site (the right boundary) of the genomic block with Eqs.62 and 63
#' @param uncollapsed_bwd_path_lbndry the uncollapsed backward path at the last site (the right boundary) of the genomic block
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @param M_sample the number of the backward paths generated
#' @return the full-collapsed backward path at the last site (the right boundary) of the genomic block returned as a matrix
collapseBwdPath_Diplo_FullRed_Block <- function(uncollapsed_bwd_path_rbndry, collapsed_info, M_sample) {
  collapsed_bwd_path_rbndry <- matrix(NA, 2, M_sample)
  if (!is.matrix(uncollapsed_bwd_path_rbndry)) {
    uncollapsed_bwd_path_rbndry <- t(as.matrix(uncollapsed_bwd_path_rbndry))
  }
  for (m in 1:M_sample) {
    collapsed_bwd_path_rbndry[1, m] <- collapsed_info$map[uncollapsed == uncollapsed_bwd_path_rbndry[1, m]]$collapsed
    collapsed_bwd_path_rbndry[2, m] <- collapsed_info$map[uncollapsed == uncollapsed_bwd_path_rbndry[2, m]]$collapsed
  }
  
  return(collapsed_bwd_path_rbndry)
}

#' Update the full-collapsed backward path at the current site of the genomic block with Eqs.64 and 65
#' @param collapsed_bwd_path_pre the full-collapsed backward path at the preceding site of the genomic block
#' @param collapsed_fwd_prob the full-collapsed forward probabilities at the current site of the genomic block 
#' @param theta the switch rate between the preceding and current sites of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param N_block the numbers of the duplicates of the haplotypes in the genomic block
#' @param M_sample the number of the backward paths generated
#' @return the full-collapsed backward path at the current site of the genomic block returned as a list
updateBwdPath_Diplo_FullRed_Block <- function(collapsed_bwd_path_pre, collapsed_fwd_prob, theta, N_panel, N_block, M_sample) {
  collapsed_bwd_path <- matrix(NA, 2, M_sample)
  for (m in 1:M_sample) {
    trans_prob <- calculateTransProb_Diplo(collapsed_bwd_path_pre[1, m], collapsed_bwd_path_pre[2, m], theta, N_panel, N_block)
    sampling_prob <- as.vector(collapsed_fwd_prob * trans_prob / sum(collapsed_fwd_prob * trans_prob))
    collapsed_bwd_path_index <- sample.int(length(N_block) * length(N_block), size = 1, replace = TRUE, prob = sampling_prob)
    collapsed_bwd_path[, m] <- if (collapsed_bwd_path_index %% length(N_block) == 0) c(length(N_block), collapsed_bwd_path_index %/% length(N_block)) else c(collapsed_bwd_path_index %% length(N_block), collapsed_bwd_path_index %/% length(N_block) + 1)
  }
  
  return(collapsed_bwd_path)
}

#' Uncollapse the full-collapsed backward path at the first site (the left boundary) of the genomic block with Eqs.66 and 67
#' @param collapsed_bwd_path_lbndry the full-collapsed backward path at the first site (the left boundary) of the genomic block
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @param M_sample the number of the backward paths generated
#' @return the uncollapsed backward path at the first site of the genomic block returned as a matrix
uncollapseBwdPath_Diplo_FullRed_Block <- function(collapsed_bwd_path_lbndry, collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, collapsed_info, M_sample) {
  uncollapsed_bwd_path_lbndry <- matrix(NA, 2, M_sample)
  for (m in 1:M_sample) {
    uncollapsed_bwd_path_lbndry_index1 <- collapsed_info$map[collapsed == collapsed_bwd_path_lbndry[1, m]]$uncollapsed
    uncollapsed_bwd_path_lbndry_index2 <- collapsed_info$map[collapsed == collapsed_bwd_path_lbndry[2, m]]$uncollapsed
    sampling_prob <- as.vector(uncollapsed_fwd_prob_lbndry[uncollapsed_bwd_path_lbndry_index1, uncollapsed_bwd_path_lbndry_index2] / collapsed_fwd_prob_lbndry[collapsed_bwd_path_lbndry[1, m], collapsed_bwd_path_lbndry[2, m]])
    uncollapsed_bwd_path_lbndry_index <- sample.int(length(uncollapsed_bwd_path_lbndry_index1) * length(uncollapsed_bwd_path_lbndry_index2), size = 1, replace = TRUE, prob = sampling_prob)
    uncollapsed_bwd_path_lbndry[1, m] <- uncollapsed_bwd_path_lbndry_index1[if (uncollapsed_bwd_path_lbndry_index %% length(uncollapsed_bwd_path_lbndry_index1) == 0) length(uncollapsed_bwd_path_lbndry_index1) else uncollapsed_bwd_path_lbndry_index %% length(uncollapsed_bwd_path_lbndry_index1)]
    uncollapsed_bwd_path_lbndry[2, m] <- uncollapsed_bwd_path_lbndry_index2[if (uncollapsed_bwd_path_lbndry_index %% length(uncollapsed_bwd_path_lbndry_index1) == 0) uncollapsed_bwd_path_lbndry_index %/% length(uncollapsed_bwd_path_lbndry_index1) else uncollapsed_bwd_path_lbndry_index %/% length(uncollapsed_bwd_path_lbndry_index1) + 1]
  }
  
  return(uncollapsed_bwd_path_lbndry)
}

#' Perform the backward procedure using the haplotype reference block with the full-collapsed state space
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param lambda the mutation rate
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_bwd_path_rbndry the uncollapsed backward path at the last site (the right boundary) of the genomic block
#' @param diplotype_rbndry the diplotype at the last site (the right boundary) of the genomic block
#' @param collapsed_fwd_prob the full-collapsed forward probabilities at the all sites of the genomic block 
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param M_sample the number of the backward paths generated
#' @return the full-collapsed backward path and the diplotype across the genomic block returned as a list
performBwdSampling_Diplo_FullRed_Block <- function(obs_block, uncollapsed_ref_block, theta_block, lambda, N_panel, uncollapsed_bwd_path_rbndry, diplotype_rbndry, collapsed_fwd_prob, uncollapsed_fwd_prob_lbndry, M_sample) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_bwd_path <- list()
  diplotype <- array(NA, dim = c(2, L_block, M_sample))
  diplotype[, L_block, ] <- diplotype_rbndry
  
  # generate the full-collapsed backward path at the last site of the genomic block
  collapsed_bwd_path[[L_block]] <- collapseBwdPath_Diplo_FullRed_Block(uncollapsed_bwd_path_rbndry, collapsed_info, M_sample)
  
  if (L_block > 1) {
    # generate the full-collapsed backward path at the following sites of the genomic block
    for (l in (L_block - 1):1) {
      collapsed_bwd_path[[l]] <- updateBwdPath_Diplo_FullRed_Block(collapsed_bwd_path[[l + 1]], collapsed_fwd_prob[[l]]$total, theta_block[l], N_panel, N_block, M_sample)
      for (m in 1:M_sample) {
        ref_haplotype1 <- collapsed_ref_block[collapsed_bwd_path[[l]][1, m], l]
        ref_haplotype2 <- collapsed_ref_block[collapsed_bwd_path[[l]][2, m], l]
        diplotype[, l, m] <- generateDiplotype_Diplo(obs_block[l], ref_haplotype1, ref_haplotype2, lambda)
      }
    }
  }
  
  # generate the uncollapsed backward path at the last site of the genomic block
  uncollapsed_bwd_path_lbndry <- uncollapseBwdPath_Diplo_FullRed_Block(collapsed_bwd_path[[1]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, collapsed_info, M_sample)
  
  return(list(collapsed_bwd_path = collapsed_bwd_path, uncollapsed_bwd_path_lbndry = uncollapsed_bwd_path_lbndry, diplotype = diplotype))
}

#########################

#' Perform the backward procedure using the haplotype reference panel with state-space reduction
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param lambda the mutation rate
#' @param block_length the length of the genomic block
#' @param collapsed_fwd_prob the full-collapsed forward probabilities across the whole genome
#' @param uncollapsed_fwd_prob_bndry the uncollapsed forward probabilities at the boundaries of the genome blocks
#' @param M_sample the number of the backward paths generated
#' @return the full-collapsed backward path across the whole genome, the uncollapsed backward path at the boundaries of the genome blocks and the diplotype across the whole genome returned as a list
performBwdProcedure_Diplo_Red_Panel <- function(obs_individual, ref_panel, theta, lambda, block_length, collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, M_sample) {
  N_panel <- nrow(ref_panel)
  L_panel <- ncol(ref_panel)
  
  # execute the preprocessing
  obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
  ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
  theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
  
  K_panel <- length(ref_panel_partition)
  
  # execute the backward procedure
  bwd_procedure_initialisation <- initialiseBwdPath_Diplo_Red_Panel(obs_individual[L_panel], ref_panel[, L_panel], lambda, N_panel, uncollapsed_fwd_prob_bndry[[K_panel + 1]], M_sample)
  uncollapsed_bwd_path_bndry <- list()
  uncollapsed_bwd_path_bndry[[K_panel + 1]] <- bwd_procedure_initialisation$bwd_path
  diplotype_block <- list()
  diplotype_block[[K_panel + 1]] <- bwd_procedure_initialisation$diplotype
  
  collapsed_bwd_path <- list()
  
  bwd_procedure_block <- performBwdSampling_Diplo_FullRed_Block(obs_individual_partition[[K_panel]], ref_panel_partition[[K_panel]], theta_partition[[K_panel]], lambda, N_panel, uncollapsed_bwd_path_bndry[[K_panel + 1]], diplotype_block[[K_panel + 1]], collapsed_fwd_prob[[K_panel]], uncollapsed_fwd_prob_bndry[[K_panel]], M_sample)
  collapsed_bwd_path[[K_panel]] <- aperm(array(unlist(bwd_procedure_block$collapsed_bwd_path), dim = c(2, M_sample, length(obs_individual_partition[[K_panel]]))), c(1, 3, 2))
  uncollapsed_bwd_path_bndry[[K_panel]] <- bwd_procedure_block$uncollapsed_bwd_path_lbndry
  diplotype_block[[K_panel]] <- bwd_procedure_block$diplotype
  
  if (K_panel > 1) {
    for (k in (K_panel - 1):1) {
      bwd_procedure_block <- performBwdSampling_Diplo_FullRed_Block(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], lambda, N_panel, uncollapsed_bwd_path_bndry[[k + 1]], diplotype_block[[k + 1]][, 1, ], collapsed_fwd_prob[[k]], uncollapsed_fwd_prob_bndry[[k]], M_sample)
      collapsed_bwd_path[[k]] <- aperm(array(unlist(bwd_procedure_block$collapsed_bwd_path), dim = c(2, M_sample, length(obs_individual_partition[[k]]))), c(1, 3, 2))
      uncollapsed_bwd_path_bndry[[k]] <- bwd_procedure_block$uncollapsed_bwd_path_lbndry
      diplotype_block[[k]] <- bwd_procedure_block$diplotype
    }
  }
  
  #
  diplotype_panel <- array(NA, dim = c(2, L_panel, M_sample))
  k_flag <- 1
  for (k in 1:K_panel) {
    L_block <- dim(diplotype_block[[k]])[2]
    diplotype_panel[, k_flag:(k_flag + L_block - 1), ] <- diplotype_block[[k]]
    k_flag <- k_flag + L_block - 1
  }
  
  return(list(collapsed_bwd_path = collapsed_bwd_path, uncollapsed_bwd_path_bndry = uncollapsed_bwd_path_bndry, diplotype_block = diplotype_block[-(K_panel + 1)], diplotype_panel = diplotype_panel))
}

###########################################################################

#' @section Section 2.2.2 Diplotype distribution

#' Construct the empirical diplotype distribution
#' @param diplo_sample the sampled diplotypes across the whole genome
#' @return the diplotypes across the whole genome and their frequencies returned as a list
constructDiploDistribution <- function(diplo_sample) {
  M_sample <- dim(diplo_sample)[3]
  
  diplotype <- list()
  frequency <- list()
  
  index_range <- 1:M_sample
  n <- 1
  
  while (length(index_range) > 0) {
    index <- NULL
    for (m in index_range) {
      if (all(diplo_sample[, , m] == diplo_sample[1:2, , index_range[1]]) || all(diplo_sample[, , m] == diplo_sample[2:1, , index_range[1]])) {
        index <- c(index, m)
      }
    }
    diplotype[[n]] <- diplo_sample[, , index_range[1]]
    frequency[[n]] <- length(index) / M_sample
    
    index_range <- index_range[!(index_range %in% index)]
    n <- n + 1
  }
  
  return(list(diplotype = diplotype, frequency = frequency))
}

##################################################

#' Perform the block-wise path sampling method
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param Ne the effective population size
#' @param d the genetic distances between successive sites
#' @param block_length the length of the genomic block
#' @param M_sample the number of the backward paths generated
performBWPS <- function(obs_individual, ref_panel, theta, lambda, block_length, M_sample) {
  
  performFwdProcedure_Diplo_Red_Panel <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
    
    performFwdProcedure_Diplo_Red_Block <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
      
      collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
      
      collapsed_ref_block <- collapsed_info$unq
      N_block <- collapsed_info$cnt
      L_block <- ncol(collapsed_ref_block)
      
      collapsed_fwd_prob <- list()
      semi_collapsed_fwd_prob <- list()
      
      # calculate the full-collapsed forward probabilities at the first site of the genomic block
      collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
      # calculate the semi-collapsed forward probabilities at the first site of the genomic block
      semi_collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_SemiRed_Block(uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info))
      
      if (min(collapsed_fwd_prob[[1]]$total) < 1e-09 || min(semi_collapsed_fwd_prob[[1]]$total) < 1e-09) {
        collapsed_fwd_prob[[1]]$total <- collapsed_fwd_prob[[1]]$total * 1e+09
        semi_collapsed_fwd_prob[[1]]$total <- semi_collapsed_fwd_prob[[1]]$total * 1e+09
      }
      
      if (L_block > 1) {
        # calculate the semi-collapsed forward probabilities at the second site of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
        collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], uncollapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
        semi_collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_SemiRed_Block(semi_collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
        
        if (min(collapsed_fwd_prob[[2]]$nr_nr) < 1e-09 || min(semi_collapsed_fwd_prob[[2]]$nr_nr) < 1e-09) {
          collapsed_fwd_prob[[2]]$nr_nr <- collapsed_fwd_prob[[2]]$nr_nr * 1e+09
          collapsed_fwd_prob[[2]]$nr_r <- collapsed_fwd_prob[[2]]$nr_r * 1e+09
          collapsed_fwd_prob[[2]]$r_nr <- collapsed_fwd_prob[[2]]$r_nr * 1e+09
          collapsed_fwd_prob[[2]]$r_r <- collapsed_fwd_prob[[2]]$r_r * 1e+09
          collapsed_fwd_prob[[2]]$total <- collapsed_fwd_prob[[2]]$total * 1e+09
          semi_collapsed_fwd_prob[[2]]$nr_nr <- semi_collapsed_fwd_prob[[2]]$nr_nr * 1e+09
          semi_collapsed_fwd_prob[[2]]$nr_r <- semi_collapsed_fwd_prob[[2]]$nr_r * 1e+09
        }
        
        if (L_block > 2) {
          # calculate the semi-collapsed forward probabilities at the following sites of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
          for (l in 3:L_block) {
            collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
            collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
            collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], uncollapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
            semi_collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_SemiRed_Block(semi_collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
            
            if (min(collapsed_fwd_prob[[l]]$nr_nr) < 1e-09 || min(semi_collapsed_fwd_prob[[l]]$nr_nr) < 1e-09) {
              collapsed_fwd_prob[[l]]$nr_nr <- collapsed_fwd_prob[[l]]$nr_nr * 1e+09
              collapsed_fwd_prob[[l]]$nr_r <- collapsed_fwd_prob[[l]]$nr_r * 1e+09
              collapsed_fwd_prob[[l]]$r_nr <- collapsed_fwd_prob[[l]]$r_nr * 1e+09
              collapsed_fwd_prob[[l]]$r_r <- collapsed_fwd_prob[[l]]$r_r * 1e+09
              collapsed_fwd_prob[[l]]$total <- collapsed_fwd_prob[[l]]$total * 1e+09
              semi_collapsed_fwd_prob[[l]]$nr_nr <- semi_collapsed_fwd_prob[[l]]$nr_nr * 1e+09
              semi_collapsed_fwd_prob[[l]]$nr_r <- semi_collapsed_fwd_prob[[l]]$nr_r * 1e+09
            }
          }
        }
      }
      
      #' calculate the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & non-recombination case and the recombination & recombination case
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
      uncollapsed_fwd_prob_rbndry_nr_nr <- uncollapsed_fwd_prob_rbndry$nr_nr
      uncollapsed_fwd_prob_rbndry_r_r <- uncollapsed_fwd_prob_rbndry$r_r
      #' calculate the uncollapsed forward probabilities at the last site of the genomic block for the non-recombination & recombination case and the recombination & non-recombination case
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_SemiRed_Block(semi_collapsed_fwd_prob[[L_block]], N_panel, collapsed_info)
      uncollapsed_fwd_prob_rbndry_nr_r <- uncollapsed_fwd_prob_rbndry$nr_r
      uncollapsed_fwd_prob_rbndry_r_nr <- uncollapsed_fwd_prob_rbndry$r_nr
      uncollapsed_fwd_prob_rbndry <- uncollapsed_fwd_prob_rbndry_nr_nr + uncollapsed_fwd_prob_rbndry_r_r + uncollapsed_fwd_prob_rbndry_nr_r + uncollapsed_fwd_prob_rbndry_r_nr
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    }
    
    N_panel <- nrow(ref_panel)
    
    # execute the preprocessing
    obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
    ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
    theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
    
    K_panel <- length(ref_panel_partition)
    
    # execute the forward procedure
    uncollapsed_fwd_prob_bndry <- list()
    uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
    
    collapsed_fwd_prob <- list()
    for (k in 1:K_panel) {
      fwd_procedure_block <- performFwdProcedure_Diplo_Red_Block(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
      collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
      uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
    }
    
    return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
  }
  
  N_panel <- nrow(ref_panel)

  emn_prob_tab <- generateEmnProbTab(lambda)
  
  fwd_proc <- performFwdProcedure_Diplo_Red_Panel(obs_individual, ref_panel, theta, block_length, emn_prob_tab)
  collapsed_fwd_prob <- fwd_proc$collapsed_fwd_prob
  uncollapsed_fwd_prob_bndry <- fwd_proc$uncollapsed_fwd_prob_bndry
  
  bwd_proc <- performBwdProcedure_Diplo_Red_Panel(obs_individual, ref_panel, theta, lambda, block_length, collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, M_sample)
  diplo_sample <- bwd_proc$diplotype_panel
  
  diplo_distribution <- constructDiploDistribution(diplo_sample)
  diplotype <- diplo_distribution$diplotype
  frequency <- unlist(diplo_distribution$frequency)
  diplotype_max_freq <- diplotype[[which(frequency == max(frequency))[1]]]
  
  return(diplotype_max_freq)
}

####################################################################################################

#' @section Section 2.3 Computational issues

#' @section Section 2.3.1 Computational complexity and optimal allocation

###########################################################################

#' @section Section 2.3.2 Improving the computational cost with approximation

#' Uncollapse the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block with Eq.80
#' @param collapsed_fwd_prob_rbndry the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the uncollapsed forward probabilities at the last site of the genomic block returned as a matrix
uncollapseFwdProb_Diplo_FullRed_Block_Approx1 <- function(collapsed_fwd_prob_rbndry, collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info) {
  uncollapsed_fwd_prob_rbndry <- matrix(NA, N_panel, N_panel)
  for (i in 1:N_panel) {
    for (j in 1:N_panel) {
      u <- collapsed_info$map[uncollapsed == i]$collapsed
      v <- collapsed_info$map[uncollapsed == j]$collapsed
      N_block_u <- collapsed_info$cnt[u]
      N_block_v <- collapsed_info$cnt[v]
      uncollapsed_fwd_prob_rbndry[i, j] <- 
        collapsed_fwd_prob_rbndry$nr_nr[u, v] * uncollapsed_fwd_prob_lbndry[i, j] / collapsed_fwd_prob_lbndry[u, v] + 
        collapsed_fwd_prob_rbndry$nr_r[u, v] * sum(uncollapsed_fwd_prob_lbndry[i, ]) / sum(collapsed_fwd_prob_lbndry[u, ]) / N_block_v + 
        collapsed_fwd_prob_rbndry$r_nr[u, v] / N_block_u * sum(uncollapsed_fwd_prob_lbndry[, j]) / sum(collapsed_fwd_prob_lbndry[, v]) + 
        collapsed_fwd_prob_rbndry$r_r[u, v] / N_block_u / N_block_v
    }
  }
  
  return(uncollapsed_fwd_prob_rbndry)
}

#' Perform the forward procedure using the haplotype reference block with state-space reduction and approximation
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the genomic block and the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block returned as a list
performFwdProcedure_Diplo_FullRed_Block_Approx1 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_fwd_prob <- list()
  
  # calculate the full-collapsed forward probabilities at the first site of the genomic block
  collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
  
  if (L_block > 1) {
    # calculate the full-collapsed forward probabilities at the second site of the genomic block
    collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
    collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
    
    if (L_block > 2) {
      # calculate the full-collapsed forward probabilities at the following sites of the genomic block
      for (l in 3:L_block) {
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
        collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
      }
    }
  }
  
  #' calculate the uncollapsed forward probabilities at the last site of the genomic block
  uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx1(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
}

#' Perform the forward procedure using the haplotype reference panel with state-space reduction and approximation
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param block_length the length of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the whole genome and the uncollapsed forward probabilities at the boundaries of the genome blocks returned as a list
performFwdProcedure_Diplo_Red_Panel_Approx1 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
  N_panel <- nrow(ref_panel)
  
  # execute the preprocessing
  obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
  ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
  theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
  
  K_panel <- length(ref_panel_partition)
  
  # execute the forward procedure
  uncollapsed_fwd_prob_bndry <- list()
  uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
  
  collapsed_fwd_prob <- list()
  for (k in 1:K_panel) {
    fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx1(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}

##################################################

#' Perform the block-wise path sampling method with approximation Eqs.80
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param Ne the effective population size
#' @param d the genetic distances between successive sites
#' @param block_length the length of the genomic block
#' @param M_sample the number of the backward paths generated
performBWPS_Approx1 <- function(obs_individual, ref_panel, theta, lambda, block_length, M_sample) {
  
  performFwdProcedure_Diplo_Red_Panel_Approx1 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
    
    performFwdProcedure_Diplo_FullRed_Block_Approx1 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
      
      collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
      
      collapsed_ref_block <- collapsed_info$unq
      N_block <- collapsed_info$cnt
      L_block <- ncol(collapsed_ref_block)
      
      collapsed_fwd_prob <- list()
      
      # calculate the full-collapsed forward probabilities at the first site of the genomic block
      collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
      
      if (min(collapsed_fwd_prob[[1]]$total) < 1e-09) {
        collapsed_fwd_prob[[1]]$total <- collapsed_fwd_prob[[1]]$total * 1e+09
      }
      
      if (L_block > 1) {
        # calculate the full-collapsed forward probabilities at the second site of the genomic block
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
        collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
        
        if (min(collapsed_fwd_prob[[2]]$nr_nr) < 1e-09) {
          collapsed_fwd_prob[[2]]$nr_nr <- collapsed_fwd_prob[[2]]$nr_nr * 1e+09
          collapsed_fwd_prob[[2]]$nr_r <- collapsed_fwd_prob[[2]]$nr_r * 1e+09
          collapsed_fwd_prob[[2]]$r_nr <- collapsed_fwd_prob[[2]]$r_nr * 1e+09
          collapsed_fwd_prob[[2]]$r_r <- collapsed_fwd_prob[[2]]$r_r * 1e+09
          collapsed_fwd_prob[[2]]$total <- collapsed_fwd_prob[[2]]$total * 1e+09
        }
        
        if (L_block > 2) {
          # calculate the semi-collapsed forward probabilities at the following sites of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
          for (l in 3:L_block) {
            collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
            collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
            
            if (min(collapsed_fwd_prob[[l]]$nr_nr) < 1e-09) {
              collapsed_fwd_prob[[l]]$nr_nr <- collapsed_fwd_prob[[l]]$nr_nr * 1e+09
              collapsed_fwd_prob[[l]]$nr_r <- collapsed_fwd_prob[[l]]$nr_r * 1e+09
              collapsed_fwd_prob[[l]]$r_nr <- collapsed_fwd_prob[[l]]$r_nr * 1e+09
              collapsed_fwd_prob[[l]]$r_r <- collapsed_fwd_prob[[l]]$r_r * 1e+09
              collapsed_fwd_prob[[l]]$total <- collapsed_fwd_prob[[l]]$total * 1e+09
            }
          }
        }
      }
      
      #' calculate the uncollapsed forward probabilities at the last site of the genomic block
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx1(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    }
    
    N_panel <- nrow(ref_panel)
    
    # execute the preprocessing
    obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
    ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
    theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
    
    K_panel <- length(ref_panel_partition)
    
    # execute the forward procedure
    uncollapsed_fwd_prob_bndry <- list()
    uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
    
    collapsed_fwd_prob <- list()
    for (k in 1:K_panel) {
      fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx1(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
      collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
      uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
    }
    
    return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
  }
  
  N_panel <- nrow(ref_panel)
  
  emn_prob_tab <- generateEmnProbTab(lambda)
  
  fwd_proc <- performFwdProcedure_Diplo_Red_Panel_Approx1(obs_individual, ref_panel, theta, block_length, emn_prob_tab)
  collapsed_fwd_prob <- fwd_proc$collapsed_fwd_prob
  uncollapsed_fwd_prob_bndry <- fwd_proc$uncollapsed_fwd_prob_bndry
  
  bwd_proc <- performBwdProcedure_Diplo_Red_Panel(obs_individual, ref_panel, theta, lambda, block_length, collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, M_sample)
  diplo_sample <- bwd_proc$diplotype_panel
  
  diplo_distribution <- constructDiploDistribution(diplo_sample)
  diplotype <- diplo_distribution$diplotype
  frequency <- unlist(diplo_distribution$frequency)
  diplotype_max_freq <- diplotype[[which(frequency == max(frequency))[1]]]
  
  return(diplotype_max_freq)
}

##################################################

#' Uncollapse the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block with Eq.84
#' @param collapsed_fwd_prob_rbndry the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the uncollapsed forward probabilities at the last site of the genomic block returned as a matrix
uncollapseFwdProb_Diplo_FullRed_Block_Approx2 <- function(collapsed_fwd_prob_rbndry, collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info) {
  uncollapsed_fwd_prob_rbndry <- matrix(NA, N_panel, N_panel)
  for (i in 1:N_panel) {
    for (j in 1:N_panel) {
      u <- collapsed_info$map[uncollapsed == i]$collapsed
      v <- collapsed_info$map[uncollapsed == j]$collapsed
      N_block_u <- collapsed_info$cnt[u]
      N_block_v <- collapsed_info$cnt[v]
      uncollapsed_fwd_prob_rbndry[i, j] <- 
        collapsed_fwd_prob_rbndry$nr_nr[u, v] * (sum(uncollapsed_fwd_prob_lbndry[i, ]) / sum(collapsed_fwd_prob_lbndry[u, ])) * (sum(uncollapsed_fwd_prob_lbndry[, j]) / sum(collapsed_fwd_prob_lbndry[, v])) + 
        collapsed_fwd_prob_rbndry$nr_r[u, v] * sum(uncollapsed_fwd_prob_lbndry[i, ]) / sum(collapsed_fwd_prob_lbndry[u, ]) / N_block_v + 
        collapsed_fwd_prob_rbndry$r_nr[u, v] / N_block_u * sum(uncollapsed_fwd_prob_lbndry[, j]) / sum(collapsed_fwd_prob_lbndry[, v]) + 
        collapsed_fwd_prob_rbndry$r_r[u, v] / N_block_u / N_block_v
    }
  }
  
  return(uncollapsed_fwd_prob_rbndry)
}

#' Perform the forward procedure using the haplotype reference block with state-space reduction and approximation
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the genomic block and the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block returned as a list
performFwdProcedure_Diplo_FullRed_Block_Approx2 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_fwd_prob <- list()
  
  # calculate the full-collapsed forward probabilities at the first site of the genomic block
  collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
  
  if (L_block > 1) {
    # calculate the full-collapsed forward probabilities at the second site of the genomic block
    collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
    collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
    
    if (L_block > 2) {
      # calculate the full-collapsed forward probabilities at the following sites of the genomic block
      for (l in 3:L_block) {
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
        collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
      }
    }
  }
  
  #' calculate the uncollapsed forward probabilities at the last site of the genomic block
  uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx2(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
}

#' Perform the forward procedure using the haplotype reference panel with state-space reduction and approximation
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param block_length the length of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the whole genome and the uncollapsed forward probabilities at the boundaries of the genome blocks returned as a list
performFwdProcedure_Diplo_Red_Panel_Approx2 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
  N_panel <- nrow(ref_panel)
  
  # execute the preprocessing
  obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
  ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
  theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
  
  K_panel <- length(ref_panel_partition)
  
  # execute the forward procedure
  uncollapsed_fwd_prob_bndry <- list()
  uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
  
  collapsed_fwd_prob <- list()
  for (k in 1:K_panel) {
    fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx2(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}

##################################################

#' Perform the block-wise path sampling method with approximation Eqs.84
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param Ne the effective population size
#' @param d the genetic distances between successive sites
#' @param block_length the length of the genomic block
#' @param M_sample the number of the backward paths generated
performBWPS_Approx2 <- function(obs_individual, ref_panel, theta, lambda, block_length, M_sample) {
  
  performFwdProcedure_Diplo_Red_Panel_Approx2 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
    
    performFwdProcedure_Diplo_FullRed_Block_Approx2 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
      
      collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
      
      collapsed_ref_block <- collapsed_info$unq
      N_block <- collapsed_info$cnt
      L_block <- ncol(collapsed_ref_block)
      
      collapsed_fwd_prob <- list()
      
      # calculate the full-collapsed forward probabilities at the first site of the genomic block
      collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block(uncollapsed_fwd_prob_lbndry, collapsed_info))
      
      if (min(collapsed_fwd_prob[[1]]$total) < 1e-09) {
        collapsed_fwd_prob[[1]]$total <- collapsed_fwd_prob[[1]]$total * 1e+09
      }
      
      if (L_block > 1) {
        # calculate the full-collapsed forward probabilities at the second site of the genomic block
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
        collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
        
        if (min(collapsed_fwd_prob[[2]]$nr_nr) < 1e-09) {
          collapsed_fwd_prob[[2]]$nr_nr <- collapsed_fwd_prob[[2]]$nr_nr * 1e+09
          collapsed_fwd_prob[[2]]$nr_r <- collapsed_fwd_prob[[2]]$nr_r * 1e+09
          collapsed_fwd_prob[[2]]$r_nr <- collapsed_fwd_prob[[2]]$r_nr * 1e+09
          collapsed_fwd_prob[[2]]$r_r <- collapsed_fwd_prob[[2]]$r_r * 1e+09
          collapsed_fwd_prob[[2]]$total <- collapsed_fwd_prob[[2]]$total * 1e+09
        }
        
        if (L_block > 2) {
          # calculate the semi-collapsed forward probabilities at the following sites of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
          for (l in 3:L_block) {
            collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
            collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
            
            if (min(collapsed_fwd_prob[[l]]$nr_nr) < 1e-09) {
              collapsed_fwd_prob[[l]]$nr_nr <- collapsed_fwd_prob[[l]]$nr_nr * 1e+09
              collapsed_fwd_prob[[l]]$nr_r <- collapsed_fwd_prob[[l]]$nr_r * 1e+09
              collapsed_fwd_prob[[l]]$r_nr <- collapsed_fwd_prob[[l]]$r_nr * 1e+09
              collapsed_fwd_prob[[l]]$r_r <- collapsed_fwd_prob[[l]]$r_r * 1e+09
              collapsed_fwd_prob[[l]]$total <- collapsed_fwd_prob[[l]]$total * 1e+09
            }
          }
        }
      }
      
      #' calculate the uncollapsed forward probabilities at the last site of the genomic block
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx2(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
      
      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    }
    
    N_panel <- nrow(ref_panel)
    
    # execute the preprocessing
    obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
    ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
    theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
    
    K_panel <- length(ref_panel_partition)
    
    # execute the forward procedure
    uncollapsed_fwd_prob_bndry <- list()
    uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
    
    collapsed_fwd_prob <- list()
    for (k in 1:K_panel) {
      fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx2(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
      collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
      uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
    }
    
    return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
  }
  
  N_panel <- nrow(ref_panel)
    
  emn_prob_tab <- generateEmnProbTab(lambda)
  
  fwd_proc <- performFwdProcedure_Diplo_Red_Panel_Approx2(obs_individual, ref_panel, theta, block_length, emn_prob_tab)
  collapsed_fwd_prob <- fwd_proc$collapsed_fwd_prob
  uncollapsed_fwd_prob_bndry <- fwd_proc$uncollapsed_fwd_prob_bndry
  
  bwd_proc <- performBwdProcedure_Diplo_Red_Panel(obs_individual, ref_panel, theta, lambda, block_length, collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, M_sample)
  diplo_sample <- bwd_proc$diplotype_panel
  
  diplo_distribution <- constructDiploDistribution(diplo_sample)
  diplotype <- diplo_distribution$diplotype
  frequency <- unlist(diplo_distribution$frequency)
  diplotype_max_freq <- diplotype[[which(frequency == max(frequency))[1]]]
  
  return(diplotype_max_freq)
}

##################################################

#' Collapse the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block with Eq.89 and 90
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the collapsed forward probabilities at the first site (the left boundary) of the genomic block returned as a matrix
collapseFwdProb_Diplo_FullRed_Block_Approx3 <- function(uncollapsed_fwd_prob_lbndry, collapsed_info) {
  N_block <- length(collapsed_info$cnt)
  
  uncollapsed_marginal_fwd_prob_lbndry <- colSums(uncollapsed_fwd_prob_lbndry)
  
  collapsed_marginal_fwd_prob_lbndry <- numeric(N_block)
  for (u in seq(N_block)) {
    collapsed_marginal_fwd_prob_lbndry[u] <- sum(uncollapsed_marginal_fwd_prob_lbndry[collapsed_info$ind[[u]]])
  }
  
  collapsed_fwd_prob_lbndry <- collapsed_marginal_fwd_prob_lbndry %*% t(collapsed_marginal_fwd_prob_lbndry)
  
  return(collapsed_fwd_prob_lbndry)
}

#' Uncollapse the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block with Eq.88
#' @param collapsed_fwd_prob_rbndry the full-collapsed forward probabilities at the last site (the right boundary) of the genomic block
#' @param collapsed_fwd_prob_lbndry the full-collapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param collapsed_info the collapsed reference block with the associated map to the uncollapsed reference block
#' @return the uncollapsed forward probabilities at the last site of the genomic block returned as a matrix
uncollapseFwdProb_Diplo_FullRed_Block_Approx3 <- function(collapsed_fwd_prob_rbndry, collapsed_fwd_prob_lbndry, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info) {
  uncollapsed_marginal_fwd_prob_rbndry <- numeric(N_panel)
  for (i in 1:N_panel) {
    u <- collapsed_info$map[uncollapsed == i]$collapsed
    N_block_u <- collapsed_info$cnt[u]
    uncollapsed_marginal_fwd_prob_rbndry[i] <- 
      sum(collapsed_fwd_prob_rbndry$nr_nr[u, ] + collapsed_fwd_prob_rbndry$nr_r[u, ]) * sum(uncollapsed_fwd_prob_lbndry[i, ]) / sum(collapsed_fwd_prob_lbndry[u, ]) + 
      sum(collapsed_fwd_prob_rbndry$r_nr[u, ] + collapsed_fwd_prob_rbndry$r_r[u, ]) / N_block_u
  }

  uncollapsed_fwd_prob_rbndry <- uncollapsed_marginal_fwd_prob_rbndry %*% t(uncollapsed_marginal_fwd_prob_rbndry)
  
  return(uncollapsed_fwd_prob_rbndry)
}

#' Perform the forward procedure using the haplotype reference block with state-space reduction and approximation
#' @param obs_block the observed block
#' @param uncollapsed_ref_block the uncollapsed reference block
#' @param theta the switch rates across the genomic block
#' @param N_panel the number of the haplotypes in the reference panel
#' @param uncollapsed_fwd_prob_lbndry the uncollapsed forward probabilities at the first site (the left boundary) of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the genomic block and the uncollapsed forward probabilities at the last site (the right boundary) of the genomic block returned as a list
performFwdProcedure_Diplo_FullRed_Block_Approx3 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
  collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
  
  collapsed_ref_block <- collapsed_info$unq
  N_block <- collapsed_info$cnt
  L_block <- ncol(collapsed_ref_block)
  
  collapsed_fwd_prob <- list()
  
  # calculate the full-collapsed forward probabilities at the first site of the genomic block
  collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block_Approx3(uncollapsed_fwd_prob_lbndry, collapsed_info))
  
  if (L_block > 1) {
    # calculate the full-collapsed forward probabilities at the second site of the genomic block
    collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
    collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
    
    if (L_block > 2) {
      # calculate the full-collapsed forward probabilities at the following sites of the genomic block
      for (l in 3:L_block) {
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
        collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
      }
    }
  }
  
  #' calculate the uncollapsed forward probabilities at the last site of the genomic block
  uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx3(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
}

#' Perform the forward procedure using the haplotype reference panel with state-space reduction and approximation
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param theta the switch rates across the whole genome
#' @param block_length the length of the genomic block
#' @param emn_prob_tab the emission probability table
#' @return the collapsed forward probabilities across the whole genome and the uncollapsed forward probabilities at the boundaries of the genome blocks returned as a list
performFwdProcedure_Diplo_Red_Panel_Approx3 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
  N_panel <- nrow(ref_panel)
  
  # execute the preprocessing
  obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
  ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
  theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
  
  K_panel <- length(ref_panel_partition)
  
  # execute the forward procedure
  uncollapsed_fwd_prob_bndry <- list()
  uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
  
  collapsed_fwd_prob <- list()
  for (k in 1:K_panel) {
    fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx3(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
    collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
    uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
  }
  
  return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
}

##################################################

#' Perform the block-wise path sampling method with approximation Eqs.88-90
#' @param obs_individual the observed individual
#' @param ref_panel the reference panel
#' @param Ne the effective population size
#' @param d the genetic distances between successive sites
#' @param block_length the length of the genomic block
#' @param M_sample the number of the backward paths generated
performBWPS_Approx3 <- function(obs_individual, ref_panel, theta, lambda, block_length, M_sample) {
  
  performFwdProcedure_Diplo_Red_Panel_Approx3 <- function(obs_individual, ref_panel, theta, block_length, emn_prob_tab) {
    
    performFwdProcedure_Diplo_FullRed_Block_Approx3 <- function(obs_block, uncollapsed_ref_block, theta_block, N_panel, uncollapsed_fwd_prob_lbndry, emn_prob_tab) {
      
      collapsed_info <- collapseRefBlock(uncollapsed_ref_block)
      
      collapsed_ref_block <- collapsed_info$unq
      N_block <- collapsed_info$cnt
      L_block <- ncol(collapsed_ref_block)
      
      collapsed_fwd_prob <- list()
      
      # calculate the full-collapsed forward probabilities at the first site of the genomic block
      collapsed_fwd_prob[[1]] <- list(total = collapseFwdProb_Diplo_FullRed_Block_Approx3(uncollapsed_fwd_prob_lbndry, collapsed_info))
      
      if (min(collapsed_fwd_prob[[1]]$total) < 1e-09) {
        collapsed_fwd_prob[[1]]$total <- collapsed_fwd_prob[[1]]$total * 1e+09
      }
      if (max(collapsed_fwd_prob[[1]]$total) > 1e+09) {
        collapsed_fwd_prob[[1]]$total <- collapsed_fwd_prob[[1]]$total * 1e-09
      }
      
      if (L_block > 1) {
        # calculate the full-collapsed forward probabilities at the second site of the genomic block
        collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[2], collapsed_ref_block[, 2], collapsed_ref_block[, 2], emn_prob_tab)
        collapsed_fwd_prob[[2]] <- initialiseFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[1]]$total, collapsed_emn_prob, theta_block[1], N_panel, N_block)
        
        if (min(collapsed_fwd_prob[[2]]$nr_nr) < 1e-09) {
          collapsed_fwd_prob[[2]]$nr_nr <- collapsed_fwd_prob[[2]]$nr_nr * 1e+09
          collapsed_fwd_prob[[2]]$nr_r <- collapsed_fwd_prob[[2]]$nr_r * 1e+09
          collapsed_fwd_prob[[2]]$r_nr <- collapsed_fwd_prob[[2]]$r_nr * 1e+09
          collapsed_fwd_prob[[2]]$r_r <- collapsed_fwd_prob[[2]]$r_r * 1e+09
          collapsed_fwd_prob[[2]]$total <- collapsed_fwd_prob[[2]]$total * 1e+09
        }
        if (max(collapsed_fwd_prob[[2]]$nr_nr) > 1e+09) {
          collapsed_fwd_prob[[2]]$nr_nr <- collapsed_fwd_prob[[2]]$nr_nr * 1e-09
          collapsed_fwd_prob[[2]]$nr_r <- collapsed_fwd_prob[[2]]$nr_r * 1e-09
          collapsed_fwd_prob[[2]]$r_nr <- collapsed_fwd_prob[[2]]$r_nr * 1e-09
          collapsed_fwd_prob[[2]]$r_r <- collapsed_fwd_prob[[2]]$r_r * 1e-09
          collapsed_fwd_prob[[2]]$total <- collapsed_fwd_prob[[2]]$total * 1e-09
        }
        
        if (L_block > 2) {
          # calculate the semi-collapsed forward probabilities at the following sites of the genomic block for the non-recombination & non-recombination case and the non-recombination & recombination case
          for (l in 3:L_block) {
            collapsed_emn_prob <- calculateEmnProb_Diplo(obs_block[l], collapsed_ref_block[, l], collapsed_ref_block[, l], emn_prob_tab)
            collapsed_fwd_prob[[l]] <- updateFwdProb_Diplo_FullRed_Block(collapsed_fwd_prob[[l - 1]], collapsed_emn_prob, theta_block[l - 1], N_panel, N_block)
            
            if (max(collapsed_fwd_prob[[l]]$nr_nr) < 1e-09) {
              collapsed_fwd_prob[[l]]$nr_nr <- collapsed_fwd_prob[[l]]$nr_nr * 1e+09
              collapsed_fwd_prob[[l]]$nr_r <- collapsed_fwd_prob[[l]]$nr_r * 1e+09
              collapsed_fwd_prob[[l]]$r_nr <- collapsed_fwd_prob[[l]]$r_nr * 1e+09
              collapsed_fwd_prob[[l]]$r_r <- collapsed_fwd_prob[[l]]$r_r * 1e+09
              collapsed_fwd_prob[[l]]$total <- collapsed_fwd_prob[[l]]$total * 1e+09
            }
            if (max(collapsed_fwd_prob[[l]]$nr_nr) > 1e+09) {
              collapsed_fwd_prob[[l]]$nr_nr <- collapsed_fwd_prob[[l]]$nr_nr * 1e-09
              collapsed_fwd_prob[[l]]$nr_r <- collapsed_fwd_prob[[l]]$nr_r * 1e-09
              collapsed_fwd_prob[[l]]$r_nr <- collapsed_fwd_prob[[l]]$r_nr * 1e-09
              collapsed_fwd_prob[[l]]$r_r <- collapsed_fwd_prob[[l]]$r_r * 1e-09
              collapsed_fwd_prob[[l]]$total <- collapsed_fwd_prob[[l]]$total * 1e-09
            }
          }
        }
      }
      
      #' calculate the uncollapsed forward probabilities at the last site of the genomic block
      uncollapsed_fwd_prob_rbndry <- uncollapseFwdProb_Diplo_FullRed_Block_Approx3(collapsed_fwd_prob[[L_block]], collapsed_fwd_prob[[1]]$total, uncollapsed_fwd_prob_lbndry, N_panel, collapsed_info)

      return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_rbndry = uncollapsed_fwd_prob_rbndry))
    }
    
    N_panel <- nrow(ref_panel)
    
    # execute the preprocessing
    obs_individual_partition <- partitionObsIndividual(obs_individual, block_length, block_index = NULL)
    ref_panel_partition <- partitionRefPanel(ref_panel, block_length, block_index = NULL)
    theta_partition <- partitionSwitchRate(theta, block_length, block_index = NULL)
    
    K_panel <- length(ref_panel_partition)
    
    # execute the forward procedure
    uncollapsed_fwd_prob_bndry <- list()
    uncollapsed_fwd_prob_bndry[[1]] <- initialiseFwdProb_Diplo_Red_Panel(obs_individual_partition[[1]][1], ref_panel_partition[[1]][, 1], N_panel, emn_prob_tab)
    
    collapsed_fwd_prob <- list()
    for (k in 1:K_panel) {
      fwd_procedure_block <- performFwdProcedure_Diplo_FullRed_Block_Approx3(obs_individual_partition[[k]], ref_panel_partition[[k]], theta_partition[[k]], N_panel, uncollapsed_fwd_prob_bndry[[k]], emn_prob_tab)
      collapsed_fwd_prob[[k]] <- fwd_procedure_block$collapsed_fwd_prob
      uncollapsed_fwd_prob_bndry[[k + 1]] <- fwd_procedure_block$uncollapsed_fwd_prob_rbndry
    }
    
    return(list(collapsed_fwd_prob = collapsed_fwd_prob, uncollapsed_fwd_prob_bndry = uncollapsed_fwd_prob_bndry))
  }
  
  N_panel <- nrow(ref_panel)

  emn_prob_tab <- generateEmnProbTab(lambda)
  
  fwd_proc <- performFwdProcedure_Diplo_Red_Panel_Approx3(obs_individual, ref_panel, theta, block_length, emn_prob_tab)
  collapsed_fwd_prob <- fwd_proc$collapsed_fwd_prob
  uncollapsed_fwd_prob_bndry <- fwd_proc$uncollapsed_fwd_prob_bndry
  
  bwd_proc <- performBwdProcedure_Diplo_Red_Panel(obs_individual, ref_panel, theta, lambda, block_length, collapsed_fwd_prob, uncollapsed_fwd_prob_bndry, M_sample)
  diplo_sample <- bwd_proc$diplotype_panel
  
  diplo_distribution <- constructDiploDistribution(diplo_sample)
  diplotype <- diplo_distribution$diplotype
  frequency <- unlist(diplo_distribution$frequency)
  diplotype_max_freq <- diplotype[[which(frequency == max(frequency))[1]]]
  
  return(diplotype_max_freq)
}

#############################################################################################################################

#' @section Section 3 Results

#' @section Section 4 Discussion

#' @section Appendix

#############################################################################################################################


