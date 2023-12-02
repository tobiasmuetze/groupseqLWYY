#################################################################################################
# Author: Tobias Muetze
# Description: 
#   This script illustrates how to use the analyze_lwyy() function
#   to calculate the correlation matrix of test statistics from different analyses
#   in a group sequential design (for a single data set). 
#   The test statistic is the Wald statistics in the LWYY model.
#   The theory is described in Mütze et al. (2020; https://doi.org/10.1177/0962280218780538)
##################################################################################################

library(here)
library(Rcpp)
library(survival)
library(tidyverse)

# Source cpp function
sourceCpp(file = here("calc_corr.cpp"))


# Load example data 
study_final <- readRDS(here("exampleData.RDS"))
tibRecEvt <- as_tibble(study_final)


# Get total sample size
nTotal <- tibRecEvt$id %>% unique %>% length

# Define calendar times for data looks (e.g., after every year)
t_analysis <- 1:5 
# Number of analyses/data looks
k <- length(t_analysis)
# Initialize variables for saving results
Bi <- matrix(0, nrow = nTotal, ncol = length(t_analysis))
A <- beta_est <- robustSE <- beta_estSurv <- robustSESurv  <- info <- numeric(k)

# Perform analysis at each look 
for(i in seq_along(t_analysis)) {
  
  t1 <- t_analysis[i]
  
  # Get the data set at a specific
  # 1. Remove all the subject recruited at/after t1
  # 2. Replace the status of events happening after t1 with 0=censored
  # 3. Replace the study stop time of events happening after t1 with t1-recruit_time
  # 4. Replace the event/stoptime of events happening after t1 with t1
  study_t1 <- tibRecEvt %>% 
    filter(start_calendar < t1) %>% 
    mutate(status = replace(status, which(stop_calendar > t1), 0),
           stop_study = replace(stop_study, which(stop_calendar > t1), t1 - recruit_time[which(stop_calendar > t1)]),
           stop_calendar = replace(stop_calendar, which(stop_calendar > t1), t1))
  
  
  mat_study_t1 <- as.matrix(study_t1)
  # Run LWYY analysis on interim data using LWYY implementation in cpp
  out_lwyy <- analyze_lwyy(x = mat_study_t1, 
                           n_id = as.integer(nTotal), 
                           n_total = as.integer(nTotal))
  Bi[,i] <- out_lwyy$id_mat[, "Bi"]
  A[i] <- out_lwyy$A
  beta_est[i] <- out_lwyy$beta_est
  robustSE[i] <- out_lwyy$robustSE
  info[i] <- out_lwyy$information
  
  # Run LWYY analysis in coxph
  survOut <-coxph(Surv(time = start_study, 
                       time2 = stop_study, 
                       event = status)~group + cluster(id), 
                  data = study_t1)
  beta_estSurv[i] <- survOut$coefficients[[1]]
  robustSESurv[i] <- sqrt(survOut$var)
}

# Compare point estimates between cpp implementation and coxph
max(abs(beta_estSurv - beta_est))
# Compare SE between cpp implementation and coxph
max(abs(robustSESurv - robustSE))



# (Robust) test statistic 
teststat <- beta_est / robustSE

# Covariance matrix under the dependent increment assumption
covar_dep <- diag(1, nrow = k)
for(i in 2:k) {
  for(j in 1:(i-1)) {
    covar_dep[i, j] <- covar_dep[j, i] <- mean(Bi[,i] * Bi[,j])  / (A[i] * A[j]) / nTotal / (robustSE[i] * robustSE[j])
  }
}

# Covariance matrix under the independent increment assumption (i.e., with (i,j) = I_i / I_j)
covarInd <- sweep(x = matrix(info, nrow = k, ncol = k), MARGIN = 2, STATS = info, FUN = '/')
# Make matrix symmetric
covarInd[lower.tri(covarInd)] <- t(covarInd)[lower.tri(covarInd)]
covarInd <- sqrt(covarInd)

# Compare covariances
covar_dep-covarInd
