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

############## Section 4. Implementing estimators for real data #####################

###Exploratory analysis of the diabetesdata###

diabetes <- as.matrix(read.csv(file = 'https://raw.githubusercontent.com/altynsai/data_bases/main/diabetes.csv')[-c(9)])
head(diabetes)

#plots

columns <- colnames(diabetes)

#corplot

M = round(cor(diabetes),2)
dev.off()
corrplot(M, method = 'number', tl.col = 'Black') # colorful number


#hist
dev.off()
par(mfrow = c(2, 4))
for (i in 1:8){
  hist(diabetes[, i], xlab = colnames(diabetes)[i], main = paste('hist of ', columns[i]),
       cex.sub = 1.25, cex.main = 1.5, cex.lab = 1.25, cex.axis = 1.25)
}

#boxplot
dev.off()
par(mfrow = c(2, 4))
for (i in 1:8){
  boxplot(diabetes[, i], xlab = colnames(diabetes)[i], main = paste('boxplot of ', columns[i]),
          cex.sub = 1.25, cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25)
}


###Comparing estimators for mean and covariance estimators###

p <- ncol(diabetes)
n_obs <- nrow(diabetes)
columns_name <- colnames(diabetes)

#Standard mean and cov estimators
center_mean <- colMeans(diabetes)
cov_mean <- cov(diabetes)

#Standard SD estimator
outl <- outlyingness(diabetes, options = list(ndir = 200*p, type = 'Shift'))
sum(outl$flagZ) #flagged as outliers
weights <- weight_calc(outl$outlyingnessX, p)
center_sd <- weights%*%diabetes/sum(weights)
cov_sd <- t(diabetes - t(matrix(rep(center_sd,n_obs), nrow=p)))%*%diag(weights)%*% (diabetes - t(matrix(rep(center_sd,n_obs), nrow=p)))/sum(weights)

#huberization of the data
diabetes_huber <- diabetes
for (col in 1:p){
  med <- median(diabetes_huber[,col])
  mad <- mad(diabetes_huber[,col], constant = 1)
  for (obs in 1:n_obs){
    diabetes_huber[obs,col] <- min(max(diabetes_huber[obs,col], -qnorm(0.975)*mad+med), qnorm(0.975)*mad+med)
  }
}

#Huberized SD estimator
huber_outl <- outlyingness(x = diabetes_huber, z = diabetes, options = list(ndir = 200*p, type = 'Shift'))
sum(huber_outl$flagZ) #flagged as outliers
weights_huber <- weight_calc(huber_outl$outlyingnessZ, p)
center_huber <- weights_huber%*%diabetes/sum(weights_huber)
cov_huber <- t(diabetes - t(matrix(rep(center_huber,n_obs), nrow=p)))%*%diag(weights_huber)%*% (diabetes - t(matrix(rep(center_huber,n_obs), nrow=p)))/sum(weights_huber)


#final table with results
estim_compare <- as.data.frame(matrix(nrow = p*2+p*(p-1)/2, ncol = 7))
estim_compare[1:p,3] <- center_mean
estim_compare[(p+1):(2*p),3] <- diag(cov_mean)
estim_compare[(2*p+1):(p*2+p*(p-1)/2),3] <- cov_mean[upper.tri(cov_mean)]

estim_compare[1:p,4] <- t(center_sd)
estim_compare[(p+1):(2*p),4] <- diag(cov_sd)
estim_compare[(2*p+1):(p*2+p*(p-1)/2),4] <- cov_sd[upper.tri(cov_sd)]

estim_compare[1:p,6] <- t(center_huber)
estim_compare[(p+1):(2*p),6] <- diag(cov_huber)
estim_compare[(2*p+1):(p*2+p*(p-1)/2),6] <- cov_huber[upper.tri(cov_huber)]

estim_compare[,5] <- (estim_compare[,4]-estim_compare[,3])/estim_compare[,3]
estim_compare[,7] <- (estim_compare[,6]-estim_compare[,3])/estim_compare[,3]


estim_compare[1:p,1] <- 'mean'
estim_compare[(p+1):(2*p),1] <- 'variance'
estim_compare[(2*p+1):(p*2+p*(p-1)/2),1] <- 'covariance'

estim_compare[1:p,2] <- columns_name
estim_compare[(p+1):(2*p),2] <- columns_name
covar_names <- matrix(paste(rep(columns_name,each=length(columns_name)),columns_name, sep=':'),ncol=8)
covar_names <- covar_names[upper.tri(covar_names)]
estim_compare[(2*p+1):(p*2+p*(p-1)/2),2] <- covar_names

estim_compare[,3:7] <- round(estim_compare[,3:7], 2)
colnames(estim_compare) <- c('','Component','Standard estimator', 'SD-estimator','Diff.','Hub. SD-estimator','Diff.')
estim_compare

write_xlsx(estim_compare,"estim_compare.xlsx")

#plotting huberized and raw data

dev.off()
par(mfrow = c(2, 2))
plot(as.data.frame(diabetes)$SkinThickness, as.data.frame(diabetes)$Glucose,
     xlim = c(0,100), ylim = c(0,200), main = 'Raw data', xlab = 'SkinThickness',
     ylab = 'Glucose')
plot(as.data.frame(diabetes_huber)$SkinThickness, as.data.frame(diabetes_huber)$Glucose,
     xlim = c(0,100), ylim = c(0,200), main = 'Huberized data', xlab = 'SkinThickness',
     ylab = FALSE)
plot(as.data.frame(diabetes)$BMI, as.data.frame(diabetes)$Pregnancies,
     xlim = c(0,70), ylim = c(0,17), xlab = 'BMI', ylab = 'Pregnancies')
plot(as.data.frame(diabetes_huber)$BMI, as.data.frame(diabetes_huber)$Pregnancies,
     xlim = c(0,70), ylim = c(0,17), xlab = 'BMI', ylab = FALSE)


