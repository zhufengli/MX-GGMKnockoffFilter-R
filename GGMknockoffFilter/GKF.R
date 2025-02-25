### GGM knockoff filter with fixed hyperparameters (Algprithm 3 in the paper)

# input: 'X': data matrix (n \times p)
#        'q': nominal FDR level
#        'offset': indicate FDR (offset=1) / mFDR (offset=0) control
#        'Z_matrix_function': function used to construct Z-statistics (measure the importance of variable/knockoff to the response)
#        '...': extra arguments of 'Z_matrix_function' besides 'Nodewise_Y_Xs_Xk_list' and 'num.cores'
#        'knockoff_method': choose from ("equi" and "sdo"), method to construct knockoffs
#        'rule': choose from ("AND" and "OR"), rule used to recover the estimated edge set from the estimated neighborhoods
#        'a': choose from (1 and 0.01), parameter used in the optimization problem
#        'W_matrix_BasedOnZ': choose from ("difference" and "max_sign"): method used to construct feature statsitics (i.e. W-statistics) based on Z-statistics
#        'num.cores': number of cores used for the parallel, the default value is 1
# output: 'E_est': estimated edge set

##########################################################

library(knockoff)        # for nodewisely constructing knockoffs
library(parallel)        # parallel when nodewisely construct knockoffs and Z-statistics

##########################################################
### Main function:
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
  E_est <- E_est_givenW_func(W_matrix, q, offset, rule, a)
  
  my_list <-list("W" = W_matrix, "T"=E_est)
  return(my_list)
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
  dif<-0
  W_matrix[i,j] <- abs(W_matrix[i,j])
  T_unsign <- cal_threshold(W_matrix, q, offset, rule, a)
  if (T_original!=T_unsign){
    #print (T_unsign)
    dif<-1
  }
  print (T_unsign)
  return(dif)
}