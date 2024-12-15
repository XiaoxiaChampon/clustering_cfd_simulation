############################################################
# Copyright 2024 Xiaoxia Champon

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
############################################################
# Purpose: Adding weekend effects in modeling the probability and simulation from real data
# Author:  Xiaoxia Champon
# Date: 11/20/2024
##############################################################
#for Moore-Penrose inverse
library(MASS)

#for getting W from 3D array X
#library(abind)

# For: profiling and visualization of profiling
#library(profvis)

# For: gam 
library(mgcv)

# For: cubicspline
#library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
#library(FADPclust)

# For: kNNdist
#library(dbscan)

# For: elbow
# devtools::install_github("ahasverus/elbow")
# library(elbow)

# For: rand.index
#library(fossil)

# For: cfda method
#library(cfda)

# For: gather method
#library(tidyverse)

# ---- For: parallelization ----
# For: foreach loop
library(foreach)

###########
library(optparse)

# Define options
option_list <- list(
  make_option(c("-j", "--jobid"), type="integer", default=123,
              help="Job Index", metavar="JOBID"),
  make_option(c("-n", "--numcpus"), type="integer", default=32,
              help="Num CPUs", metavar="NUMCPUS")
  # make_option(c("-r", "--replicas"), type="integer", default=100,
  #             help="Num Replicas", metavar="NUMREPLICAS")
)

#####need for hazel
# Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)
# 
options_jobid <- options$jobid
options_numcpus <- options$numcpus
# options_replicas <- options$replicas
# #options_subjects <- options$subjects

# options_jobid <- 1
# options_numcpus <- 9

###########
# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- TRUE
time_elapsed <- list()
if(run_parallel)
{
  print("RUNNING PARALLEL")
  
  # For: makeCluster
  library(doParallel)
  
  # For: %dorng% or registerDoRNG for reproducable parallel random number generation
  library(doRNG)
  
  if(exists("initialized_parallel") && initialized_parallel == TRUE)
  {
    parallel::stopCluster(cl = my.cluster)
  }
  # n.cores <- parallel::detectCores()
  n.cores <- options_numcpus
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), " (num cores=", n.cores, ")\n")
  initialized_parallel <- TRUE
  
  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}



#' #' Create directories 
#' if (!dir.exists("outputs")){
#'   dir.create("outputs")
#' }
#' if (!dir.exists("outputs/clustersims")){
#'   dir.create("outputs/clustersims")
#' }

#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, t * n dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n

EstimateCategFuncData_multinormial_weekend_parallel <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  
  num_indv<- ncol(W)
  timeseries_length <-  nrow(W)
  category_count <- length(unique(c(W)))
  #weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  weekend_vector <- as.factor(c(rep(c(rep(0,5*24*60/20),rep(1,2*24*60/20)),4)))[1:timeseries_length]
  Z<-NULL
  # prob<-array(0, c(num_indv, timeseries_length , category_count))
  # weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  #for (i in 1:num_indv){
  
  
  Z_P_WeekendCoef <- foreach(i = 1:num_indv , .packages = c("mgcv")) %dorng%
    
    #T_rep <- foreach(this_row = 1:5) %dorng%
    { #source("./source_code/R/data_generator.R")
      
      Z_P_WeekendCoef <- numeric(timeseries_length*(category_count+category_count-1)+category_count-1)
      print(i)
      if ( length(table(W[,i])) == category_count) {
        
        fit_binom<-gam(list(W[,i]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                            ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
        ),
        family=multinom(K=category_count-1), method = method,
        control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
        optimizer=c("outer","bfgs")) 
        #####
        #print("Dino wo")
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size)]
        #extract z
        z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        #####
        g_mul_2 <- g_design[,c(1,category_count:basis_size)+basis_size]
        coef_fit_2 <- fit_binom$coefficients[c(1,category_count:basis_size)+basis_size]
        z2 <- g_mul_2 %*% as.matrix(coef_fit_2,ncol=1)
        
        weekend_vector_coef <- fit_binom$coefficients[c(category_count-1,basis_size+category_count-1)]
        pp <- predict(fit_binom,type="response")
        p1 <- pp[,1]
        p2 <- pp[,2]
        p3 <- pp[,3]
      } else {
        if (names(table(W[,i]))[2]=="3"){
          W[,i][W[,i]==3] <- 2
          basis_size_rev <- max(min(round(min(unname(table(W[,i])[2]), sum(1-unname(table(W[,i])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[,i]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
          # print("Dino ai")
          ######################
          ####to find design matrix
          g_design <- predict(fit_binom,type = "lpmatrix")
          g_mul <- g_design[,c(1,category_count:basis_size_rev)]
          coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
          #extract z
          z2 <- g_mul %*% as.matrix(coef_fit,ncol=1)
          #####
          z1 <- rep(0,timeseries_length)
          
          weekend_vector_coef <- c(0, fit_binom$coefficients[category_count-1])
          ##########################
          p3 <- predict(fit_binom,type="response")
          p1 <- 1-p3
          p2 <- rep(0,timeseries_length)
          
        }else {
          basis_size_rev <- max(min(round(min(unname(table(W[,i])[2]), sum(1-unname(table(W[,i])[2])))/2), basis_size ), 5)
          fit_binom <- gam(W[,i]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                           family = "binomial", method = method,
                           control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                           optimizer=c("outer","bfgs"))
          #print("Dino ni")
          ######################
          ####to find design matrix
          g_design <- predict(fit_binom,type = "lpmatrix")
          g_mul <- g_design[,c(1,category_count:basis_size_rev)]
          coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
          #extract z
          z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
          z2 <- rep(0,timeseries_length)
          weekend_vector_coef <- c(fit_binom$coefficients[category_count-1],0)
          ##########################
          p2 <- predict(fit_binom,type="response")
          p1 <- 1-p2
          p3 <- rep(0,timeseries_length)
        }
      } 
      #2t*n matrix
      # Z<- cbind(Z, c(z1,z2))
      # ##find probability
      # Z_cbind=cbind(z1,z2)
      # exp_z=exp(Z_cbind)
      # denominator_p=1+exp_z[,1]+exp_z[,2]
      # p1 <- exp_z[,1]/denominator_p
      # p2 <- exp_z[,2]/denominator_p
      # p3=1/denominator_p
      #3D matrix t*n*category 
      #prob[i,,] <- cbind(p1, p2, p3)
      # 5*t +2 length
      Z_P_WeekendCoef <- c(z1,z2,p1,p2,p3, weekend_vector_coef)
      
      return(Z_P_WeekendCoef)
      
    }
  
  Z_P_WeekendCoef <- do.call(rbind, Z_P_WeekendCoef)
  #t*n
  return(list(Z1_est=t(Z_P_WeekendCoef[,1:timeseries_length]), 
              Z2_est=t(Z_P_WeekendCoef[,(1+timeseries_length):(2*timeseries_length)]),
              p1_est=t(Z_P_WeekendCoef[,(1+timeseries_length*2):(3*timeseries_length)]), 
              p2_est=t(Z_P_WeekendCoef[,(1+timeseries_length*3):(4*timeseries_length)]), 
              p3_est=t(Z_P_WeekendCoef[,(1+timeseries_length*4):(5*timeseries_length)]) ,
              weekend_vector_coef = Z_P_WeekendCoef[,(1+timeseries_length*(category_count+category_count-1)):dim(Z_P_WeekendCoef)[2]]))
}



#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, n * t dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n
EstimateCategFuncData_multinormial_weekend <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , category_count))
  weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  for (i in 1:num_indv){
    #print(i)
    cat("Individual ", i, "\n")
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
      #####
      print("Dino wo")
      ####to find design matrix
      g_design <- predict(fit_binom,type = "lpmatrix")
      g_mul <- g_design[,c(1,category_count:basis_size)]
      coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size)]
      #extract z
      z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
      #####
      g_mul_2 <- g_design[,c(1,category_count:basis_size)+basis_size]
      coef_fit_2 <- fit_binom$coefficients[c(1,category_count:basis_size)+basis_size]
      z2 <- g_mul_2 %*% as.matrix(coef_fit_2,ncol=1)
      
      weekend_vector_coef[i, ] <- fit_binom$coefficients[c(category_count-1,basis_size+category_count-1)]
      pp <- predict(fit_binom,type="response")
      p1 <- pp[,1]
      p2 <- pp[,2]
      p3 <- pp[,3]
      } else {
      if (names(table(W[i,]))[2]=="3"){
        W[i,][W[i,]==3] <- 2
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        print("Dino ai")
        ######################
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size_rev)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
        #extract z
        z2 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        #####
        z1 <- rep(0,timeseries_length)
        
        weekend_vector_coef[i, ] <- c(0, fit_binom$coefficients[category_count-1])
        p3 <- predict(fit_binom,type="response")
        p1 <- 1-p3
        p2 <- rep(0,timeseries_length)
        ##########################
      }else {
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        print("Dino ni")
        ######################
        ####to find design matrix
        g_design <- predict(fit_binom,type = "lpmatrix")
        g_mul <- g_design[,c(1,category_count:basis_size_rev)]
        coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
        #extract z
        z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
        z2 <- rep(0,timeseries_length)
        weekend_vector_coef[i, ] <- c(fit_binom$coefficients[category_count-1],0)
        ##########################
        p2 <- predict(fit_binom,type="response")
        p1 <- 1-p2
        p3 <- rep(0,timeseries_length)
      }
    } 
    #2t*n matrix
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    # exp_z=exp(Z_cbind)
    # denominator_p=1+exp_z[,1]+exp_z[,2]
    # p1 <- exp_z[,1]/denominator_p
    # p2 <- exp_z[,2]/denominator_p
    # p3=1/denominator_p
    #3D matrix t*n*category 
    prob[i,,1] <- p1
    prob[i,,2] <- p2
    prob[i,,3] <- p3
    
    #weekend_vector_coef n*2
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length + timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ,
              weekend_vector_coef = weekend_vector_coef))
}




#February 14, 2017 to March 16, 2017, Tuesday, Thursday
#first read the three values and make W
# X_i1_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i1full.RData")
# X_i2_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i2full.RData")
# X_i3_full <- load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/submission_code/X_i3full.RData")

# X_i1_full <- load("X_i1full.RData")
# X_i2_full <- load("X_i2full.RData")
# X_i3_full <- load("X_i3full.RData")
# 
# no_tweet_data <- X_i1f
# tweet_no_mention_data <- X_i2f
# tweet_mention_data <- X_i3f
# 
# rm(X_i1f)
# rm(X_i2f)
# rm(X_i3f)
# 
# array_3D <- abind(no_tweet_data, tweet_no_mention_data, tweet_mention_data, along = 3)
# W_matrix <- array_3D[,,3] * 3 +  array_3D[,,2] * 2 + array_3D[,,1]

#sample_W <- matrix(sample(c(1,2,3), 300 *10, replace = TRUE), ncol= 300, nrow= 10)
#estimate z and p
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/gam_weekends.R")
# source("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/time_track_function.R")
#source("gam_weekends.R")
source("time_track_function.R")

# indiv_one_cat <- which(apply(W_matrix , 1, function(row){length(table(row)) == 1}))
# W_matrix_final <- W_matrix[-c(indiv_one_cat),]
#load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/time_stamp.RData")
 # load("time_stamp.RData")
 # timestamps01  <- t
timeKeeperStart("Twitter all users except 14 has only 1 category, and two only has 3 one time")
#W_matrix_final[-c(3001,3374),]  #3796
# table(W_matrix_final[-c(3001,3374),][3796,])
# 
# 1    2    3 
# 1985   10    5
# Error in qr.default(if (d[1L] < d[2L]) t(z) else z) : 
#   NA/NaN/Inf in foreign function call (arg 1)

#save(W_matrix_final, timestamps01, file = "Twitter_Multinomial_W_t_Data.RData")
load("Twitter_Multinomial_W_t_Data.RData")

#Twitter_ZP__WeekendCoeff_Final_Unique_Time_P <- EstimateCategFuncData_multinormial_weekend_parallel (timestamps01, t(W_matrix_final[-c(3001,3374),][-3796,]), basis_size=25, method="ML")

Twitter_ZP__WeekendCoeff_Final_Unique_Time_P <- EstimateCategFuncData_multinormial_weekend_parallel (timestamps01, t(W_matrix_final[-c(3001,3374),][-3796,][-634,][-2267,][-3563,]), basis_size=25, method="ML")

timeKeeperNext() 
W_data_used <- W_matrix_final[-c(3001,3374),][-3796,][-634,][-2267,][-3563,]
index_exclude <- c(c(3001,3374),3796 ,634, 2267, 3563)
# save(Twitter_ZP__WeekendCoeff_Final_Unique_Time, timestamps01, 
#      W_data_used, W_matrix_final, index_exclude, file = "Twitter_ZP__WeekendCoeff_Final_Unique_Time.RData")

save(Twitter_ZP__WeekendCoeff_Final_Unique_Time_P, timestamps01, 
     W_data_used, W_matrix_final, index_exclude, file = "Twitter_ZP__WeekendCoeff_Final_Unique_Time_P.RData")
# Unravel the two variables from zp
#load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/Twitter_ZP__WeekendCoeff_Final_Unique_Time.RData")
load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/Twitter_ZP__WeekendCoeff_Final_Unique_Time_P.RData")
#t*n
Z1_est <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$Z1_est
Z2_est <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$Z2_est

p1_est <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$p1_est
p2_est <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$p2_est
p3_est <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$p3_est

#weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:2000])

weekend_vector_coef_matrix <- Twitter_ZP__WeekendCoeff_Final_Unique_Time_P$weekend_vector_coef
#save(Z1_est, Z2_est, p1_est, p2_est, p3_est, timestamps01 , index_remove, file = "Twitter_Est.RData")

cols_with_z1 <- apply(Z1_est, 2, function(x) any((x<  -2000)))
which_columns1 <- which(cols_with_z1)  # Get column indices
print(length(which_columns1))
# # 5
cols_with_z2 <- apply(Z2_est, 2, function(x) any((x< -100)))
which_columns2 <- which(cols_with_z2)  # Get column indices
print(length(which_columns2)) #7, 8
# #13
# save(Z1_est, Z2_est, p1_est, p2_est, p3_est, timestamps01 ,
#      weekend_vector_coef_matrix,
#      which_columns1, which_columns2,file = "Twitter_Est_Mul_which_columns.RData")
# #load("Twitter_Est.RData")
# load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/Twitter_Est_Mul_which_columns.RData")
# 
par(mfrow = c(1,2))
matplot(timestamps01,Z1_est[,-c(which_columns1,which_columns2)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z1", cex.lab = 1.5,cex.axis = 2,cex.main=2)
matplot(timestamps01,Z2_est[,-c(which_columns1,which_columns2)],type='l', lty=1, col="light grey",
        xlab = "Time", ylab = "Value", main = "Z2",cex.lab = 1.5,cex.axis = 2,cex.main=2)



# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est[,-c(which_columns1, which_columns2)],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1", cex.lab = 1.5,cex.axis = 2,cex.main=2)
# matplot(timestamps01,Z2_est[,-c(which_columns1, which_columns2)],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2",cex.lab = 1.5,cex.axis = 2,cex.main=2)
# 
# par(mfrow = c(1,3))
# matplot(timestamps01,p1_est[,-c(which_columns1, which_columns2)],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "p1", cex.lab = 1.5,cex.axis = 2,cex.main=2)
# matplot(timestamps01,p2_est[,-c(which_columns1, which_columns2)],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "p2", cex.lab = 1.5,cex.axis = 2,cex.main=2)
# matplot(timestamps01,p3_est[,-c(which_columns1, which_columns2)],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "p3", cex.lab = 1.5,cex.axis = 2,cex.main=2)
# #################################
# #generate new Z
extract_scores_UNIVFPCA <- function (mZ1,mZ2, tt , PVE=0.95)
{

  m<- nrow(mZ1)
  n<-ncol(mZ1)

  out1 <- refund::fpca.face(Y=t(mZ1), argvals =tt, pve = 0.99)
  out2 <- refund::fpca.face(Y=t(mZ2), argvals =tt, pve = 0.99)
  O_m1 <- matrix(0, nrow=m, ncol=out1$npc)
  O_m2  <- matrix(0, nrow=m, ncol=out2$npc)

  #construct PHI
  Phi_est0 <-  rbind(cbind(out1$efunctions*sqrt(m),O_m2 ),
                     cbind(O_m1,out2$efunctions*sqrt(m)))
  Scores0 <- cbind(out1$scores, out2$scores)
  ScoresCov_0 <- cov(Scores0 )
  oute <- eigen(ScoresCov_0)

  K<- which(cumsum( oute$values)/sum(oute$values)>=PVE)[1]
  count_iter = 0
  delta=0.01
  while (K<2 && count_iter<100) {
    count_iter = count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    K<- which(cumsum( oute$values)/sum(oute$values)>=(PVE+delta))[1]
    delta=delta+0.01
  }

  Phi_est <-  Phi_est0%*% oute$vectors[,1:K] # correct eigenfns

  mZ <- rbind(mZ1, mZ2)
  Scores_est <- t(mZ) %*%Phi_est/sqrt(m)  # they are not demeaned

  return (list(scores=Scores_est, Phi= Phi_est))
}
#simulation

eigen_score <- extract_scores_UNIVFPCA (Z1_est[,-c(which_columns1, which_columns2)],
                                        Z1_est[,-c(which_columns1, which_columns2)],
                                        timestamps01 , PVE=0.95)
# 
 save(eigen_score, timestamps01, weekend_vector_coef_matrix, file = "Twitter_eigen_weekend_Final_Dec.RData")
# 
# cumsum(c(sd(eigen_score $scores[,1]),sd(eigen_score $scores[,2])))/sum(c(sd(eigen_score $scores[,1]),sd(eigen_score $scores[,2])))
# #[1] [1] 0.8018252 1.0000000
# 3837    2
# library(MASS)
# Z_after <-  t(sqrt(length(timestamps01))*eigen_score $scores%*%ginv(eigen_score$Phi))
# 
# 
# #apply(W_matrix_final[c(which_columns1, which_columns2),],1,function(row){print(table(row))})
# 
# 
# #Z_after: 2t*n
# save(eigen_score, Z_after, which_columns1, which_columns2, Z1_est,Z2_est,
#      p1_est,p2_est,p3_est, weekend_vector_coef_matrix,
#      timestamps01, W_matrix_final, file = "Twitter_Mul_Eigen_Weekend_WhichColumns.RData" )
# 
# load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/R/acj/Twitter_Mul_Eigen_Weekend_WhichColumns.RData")
# par(mfrow = c(1,2))
# matplot(timestamps01,Z1_est[,-c(which_columns1, which_columns2)]-Z_after[1:2000,],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z1-Z1after")
# lines(timestamps01,rep(0,length(timestamps01)),col="red")
# #recover
# 
# matplot(timestamps01,Z2_est[,-c(which_columns1, which_columns2)]-Z_after[2001:4000,],type='l', lty=1, col="light grey",
#         xlab = "Time", ylab = "Value", main = "Z2-Z2after")
# lines(timestamps01,rep(0,length(timestamps01)),col="red")
# ####################################
# 
# 
# 
library(ggplot2)
library(elbow)
library(dbscan)

dbscan_cluster <- function(data=scores_K, scale_eps){

  dimz=dim(data)[2]
  if (dimz<=2){minPts=4}
  if (dimz>2){minPts=2*dimz+1}

  dist=kNNdist(data, k = minPts-1)
  #########change to max increase
  distdataelbow=data.frame(sort(dist))
  distdataelbow$index=1:(dim(data)[1])
  ipoint <- elbow(data = distdataelbow)
  epsoptimal=(ipoint$sort.dist._selected)*scale_eps

  out_dbscan <- dbscan(data, eps =epsoptimal , minPts = minPts)
  return(list(nclust=dim(table(out_dbscan$cluster)), label = out_dbscan$cluster))
}

# #3.8
# # est_dbscan
# # 0    1    2    3 
# # 29 3742   43   16
# 
# #4
# # est_dbscan
# # 0    1    2    3 
# # 28 3743   43   16
tclusterdata=data.frame(cbind(eigen_score$scores, (weekend_vector_coef_matrix[-c(which_columns1, which_columns2),])))

est_dbscan <- dbscan_cluster(data=tclusterdata,1)$label
tclusterdata$Cluster=as.factor(est_dbscan)
colnames(tclusterdata)[1:2] =c("ksi1","ksi2")
tps2 <- ggplot(tclusterdata,aes(ksi1,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("Cluster Results",'\n',"(",dim(tclusterdata)[1]," Subjects",")")) +
  xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tps2
table(est_dbscan)
# est_dbscan
# est_dbscan
# 0    1    2 
# 17 3733   87
# dim(tclusterdata)[1]
# #3863 
# 
# 
# 
# 
# ####################################
# #plot
vec0=c(which(tclusterdata$Cluster==0))
vec1=c(which(tclusterdata$Cluster==1))
vec2=c(which(tclusterdata$Cluster==2))
# vec3=c(which(tclusterdata$Cluster==3))
# 
# #2000*3683
Z1_figure <- Z1_est[,-c(which_columns1, which_columns2)]

Z2_figure <- Z2_est[,-c(which_columns1, which_columns2)]

p1_figure <- p1_est[,-c(which_columns1, which_columns2)]
p2_figure <- p2_est[,-c(which_columns1, which_columns2)]
p3_figure <- p3_est[,-c(which_columns1, which_columns2)]
# argval <- timestamps01
# 
# save(Z1_figure, Z2_figure,
#      p1_figure, p2_figure,
#      p3_figure, argval, timestamps01,argval,
#      vec0, vec1, vec2, vec3,est_dbscan,
#      tclusterdata,W_matrix_final, which_columns1,which_columns2,
#      file = "Twitter_figure_label_mul_weekend.RData")
# # min(Z1_figure[,vec3])
# # [1] 2.237043
# # > max(Z1_figure[,vec3])
# # [1] 11.3504
# load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/Twitter_new_figure_label.RData")
# ###########################################################################################
# #graph W by cluster
# 
# #W is n * t
 W_matrix_after <- W_matrix_final[-c(which_columns1, which_columns2),]
# #sample one user per cluster
# user_index <- c(sample(vec1,1),sample(vec2,1),sample(vec3,1))
# 
# apply(W_matrix_after[user_index,],1,function(x){print(table(x))})
# 
# apply(W_matrix_after[vec2,],1,function(x){print(table(x))})#37, 33
# apply(W_matrix_after[vec3,],1,function(x){print(table(x))}) #16
# 
# noquote(names(which(W_matrix_after[vec2,][33,]==3)))
# for (i in 1:25){
#   which(timestamps01==noquote(names(which(W_matrix_after[vec2,][33,]==3)))[i])
# }
# 
# 
# 

# unname(W_matrix_after[user_index[1],start_number:end_number])
# unname(W_matrix_after[vec2,start_number:end_number][33,])
# unname(W_matrix_after[user_index[3],start_number:end_number])
# 
# save(user_index, W_matrix_after, start_number, end_number, t,
#      file = "User_Cluster_W.RData")
# 
# #################graph for each user
graph_user_category <- function(index_number){
  user_data <- data.frame(
    Time = t[start_number:end_number],  # 10 days
    Category = as.factor(W_matrix_after[user_index[index_number],start_number:end_number] ) # Random categories
  )
  library(dplyr)
  # Rename categories using recode
  user_data <- user_data %>%
    mutate(Category = recode(Category, "1" = "No Tweet",
                             "2" = "Tweet Non Ref Brand",
                             "3" = "Tweet on Ref Brand"))
  user_fig <- ggplot(user_data, aes(x = Time, y = Category)) +
    geom_tile(aes(fill = Category), color = "white") +
    scale_fill_manual(values = c("No Tweet" = "red",
                                 "Tweet Non Ref Brand" = "green",
                                 "Tweet on Ref Brand" = "blue")) +
    labs(title = "", x = "", y = "") +
    theme_minimal()+
    scale_x_continuous(
      breaks = c(min(t[start_number:end_number]), max(t[start_number:end_number])),          # Specify positions for the labels
      labels = c("March 11, 2017","March 12, 2017")  # Specify corresponding labels
    ) +
    theme(legend.position = "none",
          text=element_text(size = 15))
  print( user_fig)
}



cols_with_3 <- apply(W_matrix_final[vec2,], 1, function(x) {any(x==3)})
which_columns_3 <- which(cols_with_3)  # Get column indices
print(length(which_columns_3))

#how many 3
cols_count_3 <- apply(W_matrix_final[vec2,][which_columns_3,], 1, function(x) {table(x)})
which(W_matrix_final[vec2,][which_columns_3[5],]==3)
#cols_count_3
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# 1 1997 1966 1694 1981 1793 1975 1974 1982
# 2    2   25  297   17  180   24    1    1
# 3    1    9    9    2   27    1   25   17

# 0.159480760315253 0.210013908205841 0.226240148354196 0.276773296244784 0.277236903106166 
# 245               354               389               498               499 
# 0.325452016689847 0.343532684283727  0.41029207232267 0.426518312471025 0.451553082985628 
# 603               642               786               821               875 
# 0.526657394529439 0.611033843300881 0.643949930458971 0.645340751043115 0.646267964765879 
# 1037              1219              1290              1293              1295 
# 0.710709318497914 0.711172925359295 0.726935558646268  0.79044969865554 0.827074640704682 
# 1434              1435              1469              1606              1685 
# 0.8442280945758 0.845155308298563 0.846546128882707 0.909596662030598  0.91006026889198 
# 1722              1724              1727              1866              1867 
# 0.910523875753361  0.94622160407974 
# 1868              1945

which(W_matrix_final[vec2,][which_columns_3[7],]==3)
# 0.0709318497913769  0.108020398701901  0.175243393602225  0.208159480760315  0.243857209086694 
# 54                134                279                350                427 
# 0.272137227630969  0.310152990264256  0.339360222531293  0.378303198887344  0.406119610570236 
# 488                570                633                717                777 
# 0.447844228094576  0.477515067222995  0.518776077885953  0.540101993509504  0.576263328697265 
# 867                931               1020               1066               1144 
# 0.61381548446917  0.677793231339824  0.713027352804822  0.752433936022253  0.781177561427909 
# 1225               1363               1439               1524               1586 
# 0.822438572090867  0.844691701437181  0.890125173852573  0.923968474733426  0.943439962911451 
# 1675               1723               1824               1897               1939 

start_number <- 1670
end_number <- 1730
user_index_sub <- c(sample(vec1,2),sample(vec0,2))
#123  172  613 3435  633 1634
#[1] 3249   47 3035 1119 1634  452



graph_data_set <- data.frame(W_matrix_final[c(user_index_sub[c(1,2)],
                                              user_index_sub[c(3,4)]),start_number:end_number])

# graph_data_set <- data.frame(W_matrix_final[c(2187 , 588,358, 2155)
#                                             ,start_number:end_number])

graph_data_set[5,] <- W_matrix_final[vec2,][which_columns_3[5],start_number:end_number]
graph_data_set[6,] <- W_matrix_final[vec2,][which_columns_3[7],start_number:end_number]
graph_data_set$id <- c(user_index_sub, vec2[which_columns_3][5],vec2[which_columns_3][7] )

library(reshape2)
df_long <- melt(graph_data_set, 
                id.vars = "id", 
                variable.name = "time", 
                value.name = "state")
df_long$time <- as.character(df_long$time)
df_long$time <- substr(df_long$time, 2, nchar(df_long$time))
df_long$time <- as.numeric(df_long$time)
plotData(df_long, col = c("red",  "green", "blue")) +
       labs(title = "")+
  scale_x_continuous(
         breaks = c(min(t[start_number:end_number]), max(t[start_number:end_number])),          # Specify positions for the labels
         labels = c("March 11, 2017","March 12, 2017")  # Specify corresponding labels
       ) +
  theme(legend.position = "none",
        text=element_text(size = 15))
# user_index_sub
# [1] 2187  588 2697 2699

# user_index_sub
# 2047 2439 2700 2343
# 1037  616 2700 2697
#last two user_index_sub
#[1]  395 2686 3008 2345
# graph_user_category(3)
# 
# 
# user_data <- data.frame(
#   Time = t[start_number:end_number],  # 10 days
#   Category = as.factor(W_matrix_after[user_index[3],start_number:end_number] ) # Random categories
# )
# library(dplyr)
# # Rename categories using recode
# user_data <- user_data %>%
#   mutate(Category = recode(Category, "1" = "No Tweet",
#                            "2" = "Tweet Non Ref Brand",
#                            "3" = "Tweet on Ref Brand"))
# ggplot(user_data, aes(x = Time, y = Category)) +
#   geom_tile(aes(fill = Category), color = "white") +
#   scale_fill_manual(values = c("No Tweet" = "red", 
#                                "Tweet Non Ref Brand" = "green",
#                                "Tweet on Ref Brand" = "blue")) +
#   labs(title = "", x = "", y = "") +
#   theme_minimal()+
#   scale_x_continuous(
#     breaks = c(min(t[start_number:end_number]), max(t[start_number:end_number])),          # Specify positions for the labels
#     labels = c("March 11, 2017","March 12, 2017")  # Specify corresponding labels
#   ) +
#   # theme(legend.position = "none",
#   #text=element_text(size = 15))
#   theme(legend.position = "bottom",
#         text=element_text(size = 20))+
#   guides(color = guide_legend(ncol = 2))
# 
# #one day 72 data points
# #pch=15, square, pch=16 filled circle, pch=17, triangle, pch=18 is diamond
# # Assign symbols based on y values
# 
# 
# # pch_values <- ifelse(y == 1, 18,  # Square for y <= 3
# #                      +  ifelse(y ==2 , 4,  # Filled Circle for 3 < y <= 6
# #                                +   8))  # Triangle for y > 6
# # # Plot using matplot
# # par(mfrow=c(1,3))
# # # start_number <- 10
# # # end_number <- start_number + 10
# # plot(t[start_number:end_number], t(W_matrix_after[user_index[1],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
# #         main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
# #      yaxt = "n")
# # # Add custom y-axis with only 3 values
# # # axis(2, at = c(1, 2, 3),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# # #      labels = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"))
# # axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c(""))
# # 
# # axis(1,                         # Define x-axis manually
# #      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
# #      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
# #      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c("March 11, 2017","March 12, 2017"))
# # 
# # # Add a legend
# # #legend("topright", legend = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"),
# # #       pch = pch_values, col = c(1,1,1))
# # plot(t[start_number:end_number], t(W_matrix_after[user_index[2],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
# #      main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
# #      yaxt = "n")
# # # axis(2, at = c(1, 2, 3),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# # #      labels = c("No Tweet", "Tweet Non Ref Brand", "Tweet on Ref Brand"))
# # axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c(""))
# # 
# # axis(1,                         # Define x-axis manually
# #      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
# #      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
# #      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c("March 11, 2017","March 12, 2017"))
# # 
# # plot(t[start_number:end_number], t(W_matrix_after[user_index[3],start_number:end_number]), pch = pch_values, type = "p", col = 1:3, lty = 1,
# #      main = "", xlab = "", ylab = "",xaxt = "n",ylim=c(0.75,1.25),
# #      yaxt = "n")
# # axis(2, at = c(1),  #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c(""))
# # 
# # axis(1,                         # Define x-axis manually
# #      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
# #      at = c(min(t[start_number:end_number]), max(t[start_number:end_number])),
# #      #cex.lab = 1.5,cex.axis = 2,cex.main=2,
# #      labels = c("March 11, 2017","March 12, 2017"))
# 
# #############
# #how to graph category time series
# # Example categorical time series data
# # data <- data.frame(
# #   Time = as.Date("2024-01-01") + 0:9,  # 10 days
# #   Category = sample(c("A", "B", "C"), 10, replace = TRUE)  # Random categories
# # )
# # print(data)
# # library(ggplot2)
# # 
# # ggplot(data, aes(x = Time, y = Category, group = 1)) +
# #   geom_line(aes(color = Category)) +
# #   geom_point(aes(shape = Category), size = 3) +
# #   scale_y_discrete(limits = c("A", "B", "C")) +  # Ensure consistent order
# #   labs(title = "Categorical Time Series", x = "Time", y = "Category") +
# #   theme_minimal()
# # ggplot(data, aes(x = Time, y = Category)) +
# #   geom_tile(aes(fill = Category), color = "white") +
# #   scale_fill_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) +
# #   labs(title = "Categorical Heatmap Over Time", x = "Time", y = "Category") +
# #   theme_minimal()
# # #only show xaxis at certain point with labels
# # # Plot with specific x-axis labels
# # ggplot(data, aes(x = x, y = y)) +
# #   geom_line() +
# #   scale_x_continuous(
# #     breaks = c(2, 5, 8),          # Specify positions for the labels
# #     labels = c("Two", "Five", "Eight")  # Specify corresponding labels
# #   ) +
# #   labs(title = "Custom X-axis Labels", x = "Custom Points", y = "Value") +
# #   theme_minimal()
# ###########################################################################################
# 
# 
# load("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd/Twitter_new_figure_label.RData")
# # max(Z1_figure)
# # [1] 180.3087
# # > min(Z1_figure)
# # [1] -0.3315023
# # > max(Z2_figure)
# # [1] 174.6473
# # > min(Z2_figure)
# # [1] -261.4091
logit_p <- function(p){log(p / (1 - p))}

p1_figure_final <- logit_p(p1_figure)
p2_figure_final <- logit_p(p2_figure)
p3_figure_final <- logit_p(p3_figure)
# 
# 
cols_with_p3 <- apply(p3_figure_final, 2, function(x) any((x== -Inf)))
which_columnsp3 <- which(cols_with_p3)  # Get column indices
print(length(which_columnsp3))
p3_figure_final[p3_figure_final== -Inf] = logit(0.00000001)

cols_with_p2 <- apply(p2_figure_final, 2, function(x) any((x== -Inf)))
which_columnsp2 <- which(cols_with_p2)  # Get column indices
print(length(which_columnsp2))

p2_figure_final[p2_figure_final== -Inf] = logit(0.00000001)

cols_with_p1 <- apply(p1_figure_final, 2, function(x) any((x== Inf)))
which_columnsp1 <- which(cols_with_p1)  # Get column indices
print(length(which_columnsp1))

p1_figure_final[p1_figure_final== Inf] = logit(0.99999999)
#p_3 [0, 0.3240495]
#p_2 [0,0.9582004]
#p_1 [0.0417663, 1]
# cols_with_p3 <- apply(p3_figure_final, 2, function(x) any((x== Inf)))
# which_columnsp3 <- which(cols_with_p3)  # Get column indices
# print(length(which_columnsp3))

# #[1] 3 of them  259 2545 2691
# 
z_min <- min(c(Z1_figure,Z2_figure))
z_max <- max(c(Z1_figure,Z2_figure))
p_min <- min (c (p1_figure_final, p2_figure_final, p3_figure_final[,-c(which_columnsp3)]))
p_max <- max (c (p1_figure_final, p2_figure_final, p3_figure_final[,-c(which_columnsp3)]))
# p3_figure_final[p3_figure_final==Inf] = p_max
# 
# p_min1 <- min (c (p1_figure_final[,c(vec2,vec3)], p2_figure_final[,c(vec2,vec3)]))
# p_max1 <- max (c (p1_figure_final[,c(vec2,vec3)], p2_figure_final[,c(vec2,vec3)]))
# 
# z_min1 <- min(c(Z1_figure[,c(vec2,vec3)],Z2_figure[,c(vec2,vec3)]))
# z_max1 <- max(c(Z1_figure[,c(vec2,vec3)],Z2_figure[,c(vec2,vec3)]))

# 
p_min1 <- min (c (p1_figure_final[,c(vec2,vec1)], p2_figure_final[,c(vec2,vec1)]))
p_max1 <- max (c (p1_figure_final[,c(vec2,vec1)], p2_figure_final[,c(vec2,vec1)]))

z_min1 <- min(c(Z1_figure[,c(vec2,vec1)],Z2_figure[,c(vec2,vec1)]))
z_max1 <- max(c(Z1_figure[,c(vec2,vec1)],Z2_figure[,c(vec2,vec1)]))

argval <- timestamps01
# 
save(argval,Z1_figure, Z2_figure,
     p1_figure, p2_figure,which_columnsp3,
     p3_figure, p1_figure_final, p2_figure_final,
     p3_figure_final, argval, timestamps01,
     vec0, vec1, vec2, 
     tclusterdata,z_min,
     z_max,p_min,p_max,z_min1,
     z_max1,p_min1,p_max1,
     file = "Twiiter_figure_logit_mul_final_Dec.RData")
# 
# setwd("/Users/xzhao17/Documents/GitHub/xc_cj_clustering_cfd")
# load("Twiiter_figure_logit.RData")
# #par(mfrow=c(4,5))
# t <- timestamps01
par(mfrow=c(1,5))
##vec1
meanp=colMeans(t(p1_figure_final)[vec1,])
matplot(t, t(t(p1_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=paste0("Log Odds of ", expression(hat("p")^1*(t))),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp2=colMeans(t(p2_figure_final)[vec1,])
matplot(t, t(t(p2_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=paste0("Log Odds of ", expression(hat("p")^2*(t))),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))


meanp3=colMeans(t(p3_figure_final)[vec1,])
matplot(t,t(t(p3_figure_final)[vec1,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=paste0("Log Odds of ", expression(hat("p")^3*(t))),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz=colMeans(t(Z1_figure)[vec1,])
matplot(t, t(t(Z1_figure)[vec1,]),
        type='l', lty=1, col="light grey",

        #main=expression(hat("Z")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz2=colMeans(t(Z2_figure)[vec1,])
matplot(t, t(t(Z2_figure)[vec1,]),
        type='l', lty=1, col="light grey",

        #main=expression(hat("Z")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

# #######################################
# 
# ##vec2
# argval=t
par(mfrow=c(1,5))
meanp=colMeans(t(p1_figure_final)[vec2,])
matplot(t, t(t(p1_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanp2=colMeans(t(p2_figure_final)[vec2,])
matplot(t, t(t(p2_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))


meanp3=colMeans(t(p3_figure_final)[vec2,])
matplot(t, t(t(p3_figure_final)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^3*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
meanz=colMeans(t(Z1_figure)[vec2,])
matplot(t, t(t(Z1_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("Z")^1*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

meanz2=colMeans(t(Z2_figure)[vec2,])
matplot(t, t(t(Z2_figure)[vec2,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("Z")^2*(t)),
        xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

# #######################################
# 
# 
# #vec3
# meanp=colMeans(t(p1_figure_final)[vec3,])
# matplot(t, t(t(p1_figure_final)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         #main=expression(hat("p")^1*(t)),
#         xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
# lines(argval,meanp ,
#       type='l', lty=2, lwd=2, col = "red")
# abline(h = 0, col = "blue", 
#        type='l', lty=2, lwd=2)
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# 
# meanp2=colMeans(t(p2_figure_final)[vec3,])
# matplot(t, t(t(p2_figure_final)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         #main=expression(hat("p")^2*(t)),
#         xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
# lines(argval,meanp2 ,
#       type='l', lty=2, lwd=2, col = "red")
# abline(h = 0, col = "blue", 
#        type='l', lty=2, lwd=2)
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# 
# 
# meanp3=colMeans(t(p3_figure_final)[vec3,])
# matplot(t, t(t(p3_figure)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         #main=expression(hat("p")^3*(t)),
#         xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min1,p_max1))
# lines(argval,meanp3 ,
#       type='l', lty=2, lwd=2, col = "red")
# abline(h = 0, col = "blue", 
#        type='l', lty=2, lwd=2)
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
#      cex.lab = 2,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# meanz=colMeans(t(Z1_figure)[vec3,])
# matplot(t, t(t(Z1_figure)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         #main=expression(hat("Z")^1*(t)),
#         xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
# lines(argval,meanz ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# 
# meanz2=colMeans(t(Z2_figure)[vec3,])
# matplot(t, t(t(Z2_figure)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         # main=expression(hat("Z")^2*(t)),
#         xlab="", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min1,z_max1))
# lines(argval,meanz2 ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# 
# #######################################
# 
# #vec0
par(mfrow=c(1,5))
meanp=colMeans(t(p1_figure_final)[vec0,])
matplot(t, t(t(p1_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^1*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))


lines(argval,meanp ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)


meanp2=colMeans(t(p2_figure_final)[vec0,])
matplot(t, t(t(p2_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^2*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))

lines(argval,meanp2 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)



meanp3=colMeans(t(p3_figure_final)[vec0,])
matplot(t, t(t(p3_figure_final)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("p")^3*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(p_min,p_max))
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
lines(argval,meanp3 ,
      type='l', lty=2, lwd=2, col = "red")
abline(h = 0, col = "blue",
       type='l', lty=2, lwd=2)

meanz=colMeans(t(Z1_figure)[vec0,])
matplot(t, t(t(Z1_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        #main=expression(hat("Z")^1*(t)),
        xlab="Day", ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
lines(argval,meanz ,
      type='l', lty=2, lwd=2, col = "red")


meanz2=colMeans(t(Z2_figure)[vec0,])
matplot(t, t(t(Z2_figure)[vec0,]),
        type='l', lty=1, col="light grey",
        #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
        #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))

        # cex.lab=1.5, cex.axis=2,lty=1, lwd=3,ylim=c(-12,3),main=expression(hat("Z")^2*(t)),cex.main=2          main=mtext(bquote(" ")),
        xlab="Day",
        #main=expression(hat("Z")^2*(t)),
        ylab="",cex.lab = 1.5,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(z_min,z_max))
axis(1,                         # Define x-axis manually
     #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
     #at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
     at = c(t[1],t[400],t[800],t[1200],t[1600],t[2000]),
     cex.lab = 1.5,cex.axis = 2,cex.main=2,
     labels = c("0","6","12","18","24","30"))
lines(argval,meanz2 ,
      type='l', lty=2, lwd=2, col = "red")

# 
# #######################################
# 
# 
# ###############
# #plot just p3 for all clusters
# par(mfrow=c(1,4))
# #cluster 1, cluster 2, cluster 3, noise 
# #vec1
# meanp3=colMeans(t(p3_figure)[vec1,])
# matplot(t, t(t(p3_figure)[vec1,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         
#         xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec1,]),max(t(p3_figure)[vec1,])),main="Cluster One")
# lines(argval,meanp3 ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# #vec2
# meanp3=colMeans(t(p3_figure)[vec2,])
# matplot(t, t(t(p3_figure)[vec2,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         
#         xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec2,]),max(t(p3_figure)[vec2,])),main="Cluster Two")
# lines(argval,meanp3 ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      cex.lab = 1.5,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# #vec3
# meanp3=colMeans(t(p3_figure)[vec3,])
# matplot(t, t(t(p3_figure)[vec3,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         
#         xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec3,]),max(t(p3_figure)[vec3,])),main="Cluster Three")
# lines(argval,meanp3 ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      cex.lab = 2,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30"))
# #vec0
# meanp3=colMeans(t(p3_figure)[vec0,])
# matplot(t, t(t(p3_figure)[vec0,]),
#         type='l', lty=1, col="light grey",
#         #main=mtext(bquote("Estimated Latent Cruve "*widehat('Z')[i])),
#         #xlab="Number of Datapoints", ylab=expression(widehat('Z')[i]))
#         
#         main="Noise Group",
#         xlab="Day", ylab="",cex.lab = 2,cex.axis = 2,cex.main=2,xaxt="n",ylim=c(min(t(p3_figure)[vec0,]),max(t(p3_figure)[vec0,])))
# lines(argval,meanp3 ,
#       type='l', lty=2, lwd=2, col = "red")
# axis(1,                         # Define x-axis manually
#      #at = c(0,476/2000,952/2000,1428/2000,1904/2000),
#      at = c(0,400/2000,800/2000,1200/2000,1600/2000,2000/2000),
#      cex.lab = 2,cex.axis = 2,cex.main=2,
#      labels = c("0","6","12","18","24","30")) 
# ########################################################
# #######kmeans
library(NbClust)
reskmeansall = NbClust::NbClust(data =  tclusterdata[,-5], diss = NULL,
                                distance = "euclidean", min.nc = 2, max.nc = 5,
                                method = "kmeans",index="silhouette")
tclusterdataall=data.frame(tclusterdata)
tclusterdataall$Cluster=as.factor(reskmeansall$Best.partition)

library(ggplot2)
tpsall <- ggplot(tclusterdataall,aes(ksi1 ,ksi2,colour = Cluster)) + geom_point(aes(shape=Cluster),size=3)+ggtitle(paste0("K-means Cluster Results",'\n',"(",dim(tclusterdataall)[1]," Subjects",")")) +
  #xlab(expression('Score '* widehat(xi[i1]))) + ylab(expression('Score '* widehat(xi[i2])))
  xlab('MFPCA Score1 ') + ylab('MFPCA Score2 ')+ theme(plot.title = element_text(hjust = 0.5))+#theme(text = element_text(size = 20))+
  #theme(legend.text = element_text(size = 30))
  #theme(axis.text=element_text(size = 20),axis.title = element_text(size =30))
  theme(text=element_text(size = 20))
tpsall
table(tclusterdataall$Cluster)

# 1    2 
# 4 3833 
# 
# 
# # 
# # 
# # par(mfrow = c(1,2))
# # matplot(timestamps01,Z1_est,type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z1")
# # # matlines(timestamps01,Z_after[1:672,],type='l', lty=1, col="light blue")
# # # 
# # # 
# # matplot(timestamps01,Z2_est,type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z2")
# # #recover
# # matlines(timestamps01,Z_after[673:1344,],type='l', lty=1, col="light blue")
# 
# 
# # par(mfrow = c(1,2))
# # matplot(timestamps01,Z1_est[,-c(index_remove)][,2718],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z1")
# # matplot(timestamps01,Z2_est[,-c(index_remove)][,2718],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z2")
# 
# 
# 
# # par(mfrow = c(1,2))
# # matplot(timestamps01,Z1_est[,-c(index_remove)][,1],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z1")
# # matplot(timestamps01,Z2_est[,-c(index_remove)][,1],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z2")
# # 
# # 
# # par(mfrow = c(1,2))
# # matplot(timestamps01,Z1_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z1")
# # matplot(timestamps01,Z2_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],type='l', lty=1, col="light grey",
# #         xlab = "Time", ylab = "Value", main = "Z2")
# 
# 
# 
# 
# # 
# # options(max.print = 1000000)
# # #user 3274
# # index_z2 <- c(0)
# # for (i in 1:dim(Z2_est)[2]){
# #   if (length(which(Z2_est[,i]< -250000))>=1){
# #     print(paste0("user", i,"index", which(Z2_est[,i]< -250000)[1:5]))
# #     index_z2[i] <- i
# #   }
# #   
# # }
# # 
# # #1
# # index_z2 <- na.omit( index_z2)[-1]
# # index_z2
# # 
# # #3469, 3274, 3033
# # index_less <- c(0)
# # for (i in 1:dim(Z1_est)[2]){
# #   if (length(which(Z1_est[,i]< -200))>=1){
# #     print(paste0("user", i,"index", which(Z1_est[,i]< -200)[1:5]))
# #     index_less[i] <- i
# #   }
# #   
# # }
# # 
# # #7
# # index_less <- na.omit( index_less)[-1]
# # index_less
# # 
# # index_more <- c(0)
# # for (i in 1:dim(Z1_est)[2]){
# #   if (length(which(Z1_est[,i]> 200))>=1){
# #     print(paste0("user", i,"index", which(Z1_est[,i]> 200)))
# #     index_more[i] <- i
# #   }
# # }
# # #1
# # index_more <- na.omit( index_more)[-1]
# # index_more
# # 
# # 
# # 
# # #############
# # index_remove <- c(index_z2, index_less, index_more)[-1]
# # index_remove
# # 
# # eigen_score <- extract_scores_UNIVFPCA (Z1_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715],
# #                                         Z2_est[,-c(index_remove)][,-c(1,2718)][,-c(540,2717)][,-2715], timestamps01 , PVE=0.95)
# # 
# # max(eigen_score$scores[,1])
# # min(eigen_score$scores[,1])
# # which(eigen_score$scores[,1]==max(eigen_score$scores[,1]))
# # which(eigen_score$scores[,1]==min(eigen_score$scores[,1]))
# # 
# # which(eigen_score$scores[,2]==max(eigen_score$scores[,2]))
# # 
# # plot(eigen_score$scores[,1],eigen_score$scores[,2])
# # 
# # plot(eigen_score$scores[703,1],eigen_score$scores[703,2])
# # 
# # 
# # save(eigen_score, file = "Twitter_eigen.RData")
# 
# 
# ###latent curves
# #n*2t
# #Z_after <-  t(sqrt(length(timestamps01))*eigen_score$scores%*%ginv(eigen_score$Phi))
