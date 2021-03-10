library(mvtnorm)     # generate multivariate Gaussian date
library(knockoff)        # for nodewisely constructing knockoffs
library(parallel)  
library(doParallel)
source("OtherFuncs/Graphs.R")
source("OtherFuncs/FdrPowerGraphFunc.R")

source("GGMknockoffFilter/Nodewise_Y_Xs_Xk.R")
source("GGMknockoffFilter/Z_max_lambda.R")
source("GGMknockoffFilter/Z_coef.R")
source("GGMknockoffFilter/E_est_givenW_func.R")
#source("GGMknockoffFilter/GKF.R")
GKF <- function(X, q, offset, Z_matrix_function, ..., knockoff_method, rule, a, W_matrix_BasedOnZ, num.cores=1){
  
  # Obtain Y, design matrix X, and X_knockoffs nodewisely
  Nodewise_Y_Xs_Xk_list <- Nodewise_Y_Xs_Xk(X, knockoff_method, num.cores)
  
  # construct Z matrix
  Z_matrix_list <- Z_matrix_function(Nodewise_Y_Xs_Xk_list, num.cores, ...)
  Z_ori_matrix <- Z_matrix_list[[1]]
  Z_knock_matrix <- Z_matrix_list[[2]]
  
  # in case Z_coef_function is used and lambda_quantile is a vector
  if(ncol(Z_ori_matrix)!=nrow(Z_ori_matrix)) {stop("Error: Z_matrix is not a square matrix")}
  
  # construct W matrix
  if(W_matrix_BasedOnZ=="difference"){
    W_matrix <- Z_ori_matrix - Z_knock_matrix
  }
  if(W_matrix_BasedOnZ=="max_sign"){
    W_matrix <- sign(Z_ori_matrix-Z_knock_matrix) * pmax(Z_ori_matrix,Z_knock_matrix)
  }
  
  # compute threshold and return the estimated graph
  #E_est <- E_est_givenW_func(W_matrix, q, offset, rule, a)
  
  return(W_matrix)
}

cal_threshold <- function(W_matrix, q, offset, rule, a){
  
  if(a==1) {c_a=1.93}
  if(a==0.01) {c_a=102}
  
  p <- dim(W_matrix)[1]
  W_sorted <- Sort_W_func(W_matrix)
  
  if(rule=="AND"){
    m_max <- floor( q*(p-1)/c_a - a*offset)
  }
  if(rule=="OR"){
    m_max <- floor( 1/2 * q*(p-1)/c_a - a*offset)
  }
  
  # if m_max < 0, will never find a feasible point, so thresholds=Inf
  if(m_max < 0){
    Finial_thre_vector <- rep(Inf,p)
  }
  
  # if m_max >= 0
  if(m_max >= 0){
    T_candidate_mat <- T_candidate_mat_func(W_sorted, m_max)
    
    ### Check the feasiability using each column of 'T_candidate_mat'
    Feasbile_indicate <- 0
    for (i in (m_max+1):1) {
      
      if(nrow(T_candidate_mat) ==1){    # R programming issue: for the case that 'T_candidate_mat' is a vector (i.e. m_max=0).
        thre_vector_temp <- as.vector(T_candidate_mat)
      }else{
        thre_vector_temp <- T_candidate_mat[,i]
      }
      
      constraint_vector_temp <- Constraint_vector_func(W_matrix, thre_vector_temp, rule, offset, a)
      
      # check if feasible or not
      if(rule=="OR") { 
        num.vio <- sum(constraint_vector_temp > q/(p*c_a))
      }
      if(rule=="AND") { 
        num.vio <- sum(constraint_vector_temp > 2*q/(p*c_a) )
      }
      if(num.vio==0) {Feasbile_indicate<-1; break}
    }
    
    if(Feasbile_indicate==0) {Finial_thre_vector <- rep(Inf,p)} else {Finial_thre_vector <- thre_vector_temp}
  }
  
  return(Finial_thre_vector)
}

unsigned_threshold <- function(W_matrix, T_original,  q, offset, rule, a , i, j)
{
  W_matrix[i,j] <- abs(W_matrix[i,j])
  T_unsign <- cal_threshold(W_matrix, q, offset, rule, a)
  if (T_original!=T_unsign){return (1)}
  return(0)
}


#################################################
n = 1000          # number of observations
p = 100           # number of variables
q = 0.2           # nominal FDR level

### Band graph
b <- -0.6
k= 10
c = 10
add_minEigen <- 0.5

Result <- Band_graph(p, k, b, c, add_minEigen)
Omega <- Result[[1]]
Sigma <- Result[[2]]

mu = rep(0,p)
registerDoParallel(16)
### Generate data
for(i in 1:10){
  X <- rmvnorm(n, mu, Sigma)
  
  ### Estimate the edge set using GKF with fixed hyperparameter
  W_est <- GKF(X, q, offset=1, Z_matrix_function=Z_max_lambda, alpha=1,
               knockoff_method="equi", 
               rule="AND", a=0.01,
               W_matrix_BasedOnZ="max_sign", 
               num.cores=16)
  T <- cal_threshold(W_est, q, offset=1, rule="AND", a=0.01)
  difference<-
    foreach(i=1:p,.combine = 'cbind') %:%
    foreach(j=1:p, .combine='c') %dopar%{
      unsigned_threshold(W_est, T,  q, offset=1, rule="AND", a=0.01 , i, j)
    }
  print(sum(sum(difference)))
}

