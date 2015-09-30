#####################################################################################
## Generate_mean_hetero_data.R
## Given unbalanced settings, generate unbalanced mean heterogeneity data
## Author: Meilei
#####################################################################################
library(dplyr)

Generate_mean_hetero_data = function(pi_1, pi_2, c1 = 40, c2 = 40, t = c(20, 20, 20, 40)){
  # c1 is the size of Class 1
  # c2 is the size of Class 2
  n = c1 + c2;
  n_11 = floor(c1 * pi_1); # Total number of samples in the Class 1 from Batch 1.
  n_12 = c1 - n_11; # Total number of samples in the Class 1 from Batch 2.
  n_21 = floor(c2 * pi_2) # Total number of samples in the Class 2 from Batch 1.
  n_22 = c2 - n_21 # Total number of samples in the Class 2 from Batch 2.
  y = rep(c("Class1","Class1","Class2","Class2"), c(n_11, n_12, n_21, n_22))

  M = matrix(nrow = n, ncol = sum(t)) # In the data matrix M, samples are in the rows and feasures are in the columns. 
  colnames(M) = paste0('X', 1:sum(t))
  
  if(t[1] > 0){
    for(i in 1 : t[1]){
      M[,i] = rep(c(2, -2), c(40, 40)) + 
        rnorm(n, 0 ,1)
    }
  }
  
  if(t[2] > 0){
    for(i in (t[1] + 1) : sum(t[1:2])){
      M[,i] = rep(c(2, -2), c(40, 40)) + 
        rep(c(2,-2, 2,-2), c(n_11, n_12, n_21, n_22)) + 
        rnorm(n, 0, 1)
    }
  }
  
  if(t[3] > 0){
    for(i in (sum(t[1:2]) + 1) : sum(t[1:3])){
      M[,i] = rep(c(2,-2, 2,-2), c(n_11, n_12, n_21, n_22)) + 
        rnorm(n, 0 ,1)
    }
  }
  
  if(t[4] > 0){
    for(i in (sum(t[1:3]) + 1) : sum(t[1:4])){
      M[,i] = rnorm(n, 0 ,1)
    }
  }
  
  batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
  data.frame(y, M, batch, pi_1, pi_2, c1, c2)
}