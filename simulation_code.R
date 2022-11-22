# Robust Statistics Spring 2022 KU Leuven

# Stahel-Donoho Estimator, with(out) huberized outlyingness

# importing libraries
library(MASS)
library(dplyr)
library(rrcov)
library(mrfDepth)
library(writexl)
library(optimbase)
library(corrplot)
set.seed(1996)

rm(list=ls()) 


#################### Defining functions for further analysis ##########################
#Data simulating function
simulated_data_generator <- function(p, r, n_obs){
  # Mean centered at 0
  mu <- rep(0, p)
  # Variances of 1
  diag <- rep(1, p)
  # Setting up off-diagonal elements to be r_squared
  offdiag = rep(r, p*(p-1)/2)
  # Generating Covariance matrix from diag and off-diag elements
  m <- matrix(NA, ncol = p, nrow = p)
  diag(m) <- diag
  m[lower.tri(m)] <- offdiag
  m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
  
  # Generating data from covariance matrix above
  simulated_data <- mvrnorm(n = n_obs, mu = mu, Sigma = m)
  return(simulated_data)
}


# Inserting outliers
incorp_outs <- function(data, d, epsilon, k, sdev_contam){
  for (i in 1:d){
    numb_out <- round(nrow(data)*epsilon) 
    index_out <- sample(1:nrow(data), numb_out, replace=FALSE)
    data[index_out,i] <- rnorm(numb_out, mean=k/sqrt(d), sd=sdev_contam)
  }
  return(data)
}

# Function to calculate weights from outlyingness
weight_calc <- function(outl, p){
  c <- min(sqrt(qchisq(0.5, p)), 4)
  weights <- rep(NA, length(outl))
  for (x in 1:length(outl)){
    if (outl[x] <= c){
      weights[x] <-1
      } 
    else{
      weights[x] <- (c/outl[x])^2
    }
  }
  return (weights)
}

# Function to calculate the mse for both estimators

SD_estimators <- function(p, r, n_obs, epsilon, d, k,sdev_contam){
  mse_center <- matrix(nrow = 500, ncol = p)
  mse_diag <- matrix(nrow = 500, ncol = p)
  mse_offdiag <- matrix(nrow = 500, ncol = p*(p-1)/2)
  
  mse_center_huber <- matrix(nrow = 500, ncol = p)
  mse_diag_huber <- matrix(nrow = 500, ncol = p)
  mse_offdiag_huber <- matrix(nrow = 500, ncol = p*(p-1)/2)
  
  
  for (i in 1:500){
    data_sim <- simulated_data_generator(p, r, n_obs)
    data_sim_w_out <- incorp_outs(data_sim, d, epsilon, k, sdev_contam)
    
    #Classic Stahel-Donoho estimator 
    
    # sd_est <- CovSde(data_sim_w_out, nsamp = 200*p)
    # sd_cov <- sd_est$cov
    # mse_center[i,] <-sd_est$center^2
    # mse_diag[i,] <- (diag(sd_cov)-rep(1,p))^2
    # mse_offdiag[i,] <-(c(sd_cov[upper.tri(sd_cov)])-rep(r, p*(p-1)/2))^2
    
    outl <- outlyingness(data_sim_w_out, options = list(ndir = 200*p, type = 'Shift'))
    weights <- weight_calc(outl$outlyingnessX, p)
    center <- weights%*%data_sim_w_out/sum(weights)
    mse_center[i,] <-center^2
    cov <- t(data_sim_w_out - t(matrix(rep(center,n_obs), nrow=p)))%*%diag(weights)%*% (data_sim_w_out - t(matrix(rep(center,n_obs), nrow=p)))/sum(weights)
    mse_diag[i,] <- (diag(cov)-rep(1,p))^2
    mse_offdiag[i,] <-(c(cov[upper.tri(cov)])-rep(r,p*(p-1)/2))^2
    
    
    #Huberized data set
    data_huber <- data_sim_w_out
    for (col in 1:p){
      med <- median(data_huber[,col])
      mad <- mad(data_huber[,col], constant = 1)
      for (obs in 1:n_obs){
        data_huber[obs,col] <- min(max(data_huber[obs,col], -qnorm(0.975)*mad+med), qnorm(0.975)*mad+med)
      }
    }
    
    #Huberized SD estimator
    huber_outl <- outlyingness(x = data_huber, z = data_sim_w_out, options = list(ndir = 200*p, type = 'Shift'))
    weights <- weight_calc(huber_outl$outlyingnessZ, p)
    center_huber <- weights%*%data_sim_w_out/sum(weights)
    mse_center_huber[i,] <- (center_huber)^2
    cov_huber <- t(data_sim_w_out - t(matrix(rep(center_huber,n_obs), nrow=p)))%*%diag(weights)%*% (data_sim_w_out - t(matrix(rep(center_huber,n_obs), nrow=p)))/sum(weights)
    mse_diag_huber[i,] <- (diag(cov_huber)-rep(1,p))^2
    mse_offdiag_huber[i,] <-(c(cov_huber[upper.tri(cov_huber)])-rep(r,p*(p-1)/2))^2
  }
  
  
  #columns with contaminated components
  contam <- zeros(p,p)
  contam[1:d,1:d] <-1
  contam_2 <- as.logical(contam[upper.tri(contam)])
  
  contam <- zeros(p,p)
  contam[1:d,(d+1):p] <-1
  contam_1 <- as.logical(contam[upper.tri(contam)])
  
  #assemblying table for boxplots
  tables_for_boxplots <- cbind(rowMeans(sqrt(mse_center[,1:d])), rowMeans(sqrt(mse_center_huber[,1:d])),  rowMeans(sqrt(mse_center[,(d+1):p])), rowMeans(sqrt(mse_center_huber[,(d+1):p])))
  
  #Calculating the ratio of two estimators
  results <- c (
    mean(mse_center_huber)/mean(mse_center),
    mean(mse_center_huber[,1:d])/mean(mse_center[,1:d]),
    mean(mse_diag_huber)/mean(mse_diag),
    mean(mse_diag_huber[,1:d])/mean(mse_diag[,1:d]),
    mean(mse_offdiag_huber)/mean(mse_offdiag),
    mean(mse_offdiag_huber[,contam_1])/mean(mse_offdiag[,contam_1]),
    mean(mse_offdiag_huber[,contam_2])/mean(mse_offdiag[,contam_2])
  )
  output <- list(results, tables_for_boxplots)
  return(output)
  
  
}
###################### Section 3. Simulation data #######################
### Building the Table 1 for 5 dimensional data with epsilon 0.1 and 0.2 ### 


results <- as.data.frame(matrix(nrow = 14, ncol = 7))
table_for_boxplots <- as.data.frame(matrix(nrow = 500, ncol = 8))

output <- SD_estimators(p=5, r=0, n_obs=50, epsilon = 0.1, d=2, k=6, sdev_contam=0.1)
table_for_boxplots[,1] <- output[2]
results[1:7,4] <- output[1]

output <- SD_estimators(p=5, r=0, n_obs=50, epsilon = 0.2, d=2, k=6, sdev_contam=0.1)
table_for_boxplots[,2] <- output[2]
results[8:14,4] <- output[1]

output <- SD_estimators(p=5, r=0, n_obs=50, epsilon = 0.1, d=2, k=64, sdev_contam=0.1)
table_for_boxplots[,3] <- output[2]
results[1:7,5] <- output[1]

output <- SD_estimators(p=5, r=0, n_obs=50, epsilon = 0.2, d=2, k=64, sdev_contam=0.1)
table_for_boxplots[,4] <- output[2]
results[8:14,5] <- output[1]

output <- SD_estimators(p=5, r=0.9, n_obs=50, epsilon = 0.1, d=2, k=6, sdev_contam=0.1)
table_for_boxplots[,5] <- output[2]
results[1:7,6] <- output[1]

output <- SD_estimators(p=5, r=0.9, n_obs=50, epsilon = 0.2, d=2, k=6, sdev_contam=0.1)
table_for_boxplots[,6] <- output[2]
results[8:14,6] <- output[1]

output <- SD_estimators(p=5, r=0.9, n_obs=50, epsilon = 0.1, d=2, k=64, sdev_contam=0.1)
table_for_boxplots[,7] <- output[2]
results[1:7,7] <- output[1]

output <- SD_estimators(p=5, r=0.9, n_obs=50, epsilon = 0.2, d=2, k=64, sdev_contam=0.1)
table_for_boxplots[,8] <- output[2]
results[8:14,7] <- output[1]

results[c(1:2,8:9),1] <- 'Center'
results[c(3:4,10:11),1] <- 'Diag'
results[c(5:7,12:14),1] <- 'Offdiag'

results[1:7,2] <- 0.1
results[8:14,2] <- 0.2

results[c(1,3,5,8,10,12),3] <- 'All'
results[c(2,4,9,11),3] <- 'Cont'
results[c(6,13),3] <- '1 cont'
results[c(7,14),3] <- '2 cont'

results[,4:7] <- round(results[,4:7], digits = 2) 

colnames(results) <- c('','Epsilon', 'Components', 'r=0,k=6','r=0,k=64', 'r=0.9,k=6','r=0.9,k=64')

(results)

write_xlsx(results,"results.xlsx")

#plotting boxplots
par(mfrow=c(2,4))
boxplot(table_for_boxplots[,1], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=6, r=0')
boxplot(table_for_boxplots[,3], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=64, r=0')
boxplot(table_for_boxplots[,5], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=6, r=0.9')
boxplot(table_for_boxplots[,7], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=64, r=0.9')


boxplot(table_for_boxplots[,2], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=6, r=0')
boxplot(table_for_boxplots[,4], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=64, r=0')
boxplot(table_for_boxplots[,6], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=6, r=0.9')
boxplot(table_for_boxplots[,8], names = c('SD-C','HSD-C','SD-NC','HSD-NC'), las=2, main = 'k=64, r=0.9')

### Building the Table 2 for 7 dimensional data with epsilon 0.2 ###

results_table2 <- as.data.frame(matrix(nrow = 7, ncol = 7))

results_table2[1:7,7] <- SD_estimators(p=7, r=0.9, n_obs=100, epsilon = 0.2, d=5, k=64, sdev_contam=0.1)[1]
results_table2[1:7,6] <- SD_estimators(p=7, r=0.9, n_obs=100, epsilon = 0.2, d=5, k=6, sdev_contam=0.1)[1]
results_table2[1:7,5] <- SD_estimators(p=7, r=0, n_obs=100, epsilon = 0.2, d=5, k=64, sdev_contam=0.1)[1]
results_table2[1:7,4] <- SD_estimators(p=7, r=0, n_obs=100, epsilon = 0.2, d=5, k=6, sdev_contam=0.1)[1]

results_table2[c(1:2),1] <- 'Center'
results_table2[c(3:4),1] <- 'Diag'
results_table2[c(5:7),1] <- 'Offdiag'
results_table2[1:7,2] <- 0.2

results_table2[c(1,3,5),3] <- 'All'
results_table2[c(2,4), 3] <- 'Cont'
results_table2[c(6),3] <- '1 cont'
results_table2[c(7),3] <- '2 cont'
results_table2[,4] <- as.numeric(results_table2[,4])
results_table2[,4:7] <- round(results_table2[,4:7], digits = 2) 


colnames(results_table2) <- c('','Epsilon', 'Components', 'r=0,k=6','r=0,k=64', 'r=0.9,k=6','r=0.9,k=64')
results_table2
write_xlsx(results_table2,"results_table2.xlsx")


