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
# Purpose: Adding weekend effects in modeling the probability
# Author:  Xiaoxia Champon
# Date: 11/06/2024
##############################################################

#' Runs gam on the given data and returns results
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param in_x binary series
#' @param family_choice "binomial" or "probit"
#' @return list: fit values and linear predictors both with length of time_series length
RunGam_Day <- function(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)
{
  basis_size_rev <- max(min(round(min(sum(in_x), sum(1-in_x))/2), basis_size ), 5)
  
  if(family_choice=="binomial")
  {
    family_object <- "binomial"
  }else if(family_choice=="probit")
  {
    family_object <- binomial(link="probit")
  }
  
  fit_binom <- gam(in_x~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                   family=family_object, method = method,
                   control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                   optimizer=c("outer","bfgs"))
  prob <- fit_binom$fitted.values
  
  if(family_choice=="binomial")
  {
    prob_linpred <- fit_binom$linear.predictors
  }else if(family_choice=="probit"){
    z_probit <- fit_binom$linear.predictors
    
    prob_probit <- pnorm(z_probit)
    zprobit_reciprocal <- 1/ prob_probit
    prob_linpred <- -log(zprobit_reciprocal-1)
  }
  
  return(list(prob=prob, linpred=prob_linpred))
}

#testing code
# library(mgcv)
# #RunGam_Day <- function(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)
# timestamps01 <- seq(0,1, length = 1000)
# weekend_vector <- c(rep(c(0,0,0,0,0,1,1),200))[1:1000]
# in_x <- rbinom(1000,1,0.5)
# family_choice <- "binomial"
# basis_size <- 25
# method <- "ML"
# test_data <- RunGam_Day(timestamps01, weekend_vector, in_x, family_choice, basis_size, method)

#' Fit multinomial from W categorical data to extra smooth latent curves
#' Fits a binomial to describe the given in_x
#' @param timestamps01 current time value
#' @param W : 2D categorical data matrix, n * t dimension, n is the number of subjects, t is the time
#' @return list: fit values and linear predictors both with length of time_series length, Z is t*n, p: t*n
EstimateCategFuncData_multinormial_ind <- function(timestamps01, W, basis_size=25, method="ML")
{
  
  # num_indv<- ncol(W)
  # timeseries_length <-nrow(W)
  
  num_indv<- nrow(W)
  timeseries_length <-ncol(W)
  category_count <- length(unique(c(W)))
  weekend_vector <- as.factor(c(rep(c(rep(0,480),rep(1,192)),4))[1:timeseries_length])
  
  Z<-NULL
  prob<-array(0, c(num_indv, timeseries_length , category_count))
  weekend_vector_coef <- matrix(0, num_indv, category_count-1)
  for (i in 1:num_indv){
    print(i)
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
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
    } else {
      
      W[i,][W[i,]==3] <- 2
      
      
      basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
      fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                       family = "binomial", method = method,
                       control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                       optimizer=c("outer","bfgs"))
      
      ######################
      ####to find design matrix
      g_design <- predict(fit_binom,type = "lpmatrix")
      g_mul <- g_design[,c(1,category_count:basis_size_rev)]
      coef_fit <- fit_binom$coefficients[c(1,category_count:basis_size_rev)]
      #extract z
      z1 <- g_mul %*% as.matrix(coef_fit,ncol=1)
      #####
      # g_mul_2 <- g_design[,c(1,category_count:basis_size)+basis_size]
      # coef_fit_2 <- fit_binom$coefficients[c(1,category_count:basis_size)+basis_size]
      #z2 <- g_mul_2 %*% as.matrix(coef_fit_2,ncol=1)
      z2 <- rep(0,timeseries_length)
      
      weekend_vector_coef[i, ] <- c(fit_binom$coefficients[category_count-1],0)
      ##########################
      
      
    } 
    
    
    # z1<- fit_binom$linear.predictors[,1]
    # z2<- fit_binom$linear.predictors[,2]
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    exp_z=exp(Z_cbind)
    denominator_p=1+exp_z[,1]+exp_z[,2]
    p1 <- exp_z[,1]/denominator_p
    p2 <- exp_z[,2]/denominator_p
    p3=1/denominator_p
    prob[i,,] <- cbind(p1, p2, p3)
    
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length + timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ,
              weekend_vector_coef = weekend_vector_coef))
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
    print(i)
    if ( length(table(W[i,])) == category_count) {
      
      fit_binom<-gam(list(W[i,]-1~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector,
                          ~s(timestamps01,bs = "cc", m=2, k = basis_size) + weekend_vector
      ),
      family=multinom(K=category_count-1), method = method,
      control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
      optimizer=c("outer","bfgs")) 
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
    } else {
      if (names(table(W[i,]))[2]=="3"){
        W[i,][W[i,]==3] <- 2
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        
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
        ##########################
      }else {
        basis_size_rev <- max(min(round(min(unname(table(W[i,])[2]), sum(1-unname(table(W[i,])[2])))/2), basis_size ), 5)
        fit_binom <- gam(W[i,]-1~s(timestamps01, bs = "cc", m=2, k = basis_size_rev) + weekend_vector,
                         family = "binomial", method = method,
                         control=list(maxit = 500,mgcv.tol=1e-4,epsilon = 1e-04),
                         optimizer=c("outer","bfgs"))
        
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
      }
    } 
    #2t*n matrix
    Z<- cbind(Z, c(z1,z2))
    ##find probability
    Z_cbind=cbind(z1,z2)
    exp_z=exp(Z_cbind)
    denominator_p=1+exp_z[,1]+exp_z[,2]
    p1 <- exp_z[,1]/denominator_p
    p2 <- exp_z[,2]/denominator_p
    p3=1/denominator_p
    #3D matrix t*n*category 
    prob[i,,] <- cbind(p1, p2, p3)
    
    #weekend_vector_coef n*2
  }
  
  return(list(Z1_est=Z[1:timeseries_length ,], Z2_est=Z[1:timeseries_length + timeseries_length ,], 
              p1_est=t(prob[,,1]), p2_est=t(prob[,,2]), p3_est=t(prob[,,3]) ,
              weekend_vector_coef = weekend_vector_coef))
}

######test
#EstimateCategFuncData_multinormial <- function(timestamps01, W, basis_size=25, method="ML")
timestamps01 <- timestamps01
test_mul <- EstimateCategFuncData_multinormial_weekend (timestamps01, W_matrix_final, basis_size=25, method="ML")
