library(mvtnorm)     # generate multivariate Gaussian date
library(doParallel)
source("OtherFuncs/Graphs.R")
source("OtherFuncs/FdrPowerGraphFunc.R")

source("GGMknockoffFilter/Nodewise_Y_Xs_Xk.R")
source("GGMknockoffFilter/Z_max_lambda.R")
source("GGMknockoffFilter/Z_coef.R")
source("GGMknockoffFilter/E_est_givenW_func.R")
source("GGMknockoffFilter/GKF.R")

#################################################
n = 1200          # number of observations
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
registerDoParallel(15)

n_sim <- 100
table <- rep(0, n_sim)
### Generate data
for(i_cycle in 1:n_sim){
X <- rmvnorm(n, mu, Sigma)

### Estimate the edge set using GKF with fixed hyperparameter
my_list <- GKF(X, q, offset=1, Z_matrix_function=Z_max_lambda, alpha=1,
             knockoff_method="equi", 
             rule="AND", a=0.01,
             W_matrix_BasedOnZ="max_sign", 
             num.cores=15)
difference<-foreach(i=1:p,.combine = 'cbind') %:% foreach(j=1:p, .combine='c') %dopar% {
      unsigned_threshold(my_list$W, my_list$T,  q, offset=1, rule="AND", a=0.01 , i, j)
  }

table[i_cycle] <- sum(difference)
}


W_temp <- my_list$W
W_temp[1,71] <- abs(my_list$W[1,71])
cal_threshold(W_temp, q, offset=1, rule="AND", a=0.01)
