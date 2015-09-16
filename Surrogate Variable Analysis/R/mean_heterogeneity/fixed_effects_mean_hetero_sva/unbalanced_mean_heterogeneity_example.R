############################################################################################################
## fixed_effects_mean_heterogeneity_example.R
## This is 2-d toy example. There are two classes and each class is a Gaussian mixture. The Gaussian
## mixture will be consider as heterogeneity. SVA will be applied to alliviate this problem. SVA is expected
## to be good at regression rather than classification.
## We consider the simulation model as a two-way ANOVA setting.
##
## Author: Meilei
###########################################################################################################

library(dplyr)
library(ggplot2)
library(reshape2)


# setting of simulation ---------------------------------------------------

# Each class have 40 samples and each sample has 100 features.
# First 40 samples come from Class 1 and last 40 samples come from Class 2.
# pi_1 of Class 1 come from Batch 1 and pi_2 of Class 2 from Batch 1. 

# Primary variable (class label) effect: For Class 1, c1 = 2; For Class 2, c2 = -2.
# Batch effect: For Batch 1, b1 = 2; For Batch 2, b2 = -2.
# Random Noise: N(0, 1)

n = 80 # number of samples
t1 = 20 # Feature 1 - t1 is only affected by primary variable + random noise 
t2 = 40 # Feature t1 + 1 - t2 is affected by both primary variable and batch + random noise
t3 = 60 # Feature t2 + 1 - t3 is only affected by batch + random noise
t4 = 100 # Feature t3 + 1 - t4 is only random noise
fe = factor(rep(c('y','y & batch','batch','noise'), c(t1, t2-t1, t3-t2, t4-t3)), 
            levels = c('y','y & batch','batch','noise')) # lable of feature 
# Balanced Case -----------------------------------------------------------
pi_1 = 0.5; pi_2 = 0.5;

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)
if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

y = c(rep('Class1', 40), rep('Class2', 40))
batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 


train0 = data.frame(y, M, batch, pi_1, pi_2)

# Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()



# Unbalanced Case 1 -------------------------------------------------------

# In this case: pi_1 = pi_2, pi_1 + pi_2 < 1
# The propotion of each class in the Batch 1 are the same.
# The sample numbers of two batches are different.

## Case 1.1
pi_1 = 0.4; pi_2 = 0.4
  
n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train1.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 1.2
pi_1 = 0.1; pi_2 = 0.1

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train1.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()


# Unbalanced Case 2 -------------------------------------------------------

# pi_1 > pi_2, pi_1 + pi_2 = 1.
# The propotion of each class in the Batch 1 are different, 
# The sample numbers of two batches are the same.

## Case 2.1
pi_1 = 0.6; pi_2 = 0.4

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train2.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 2.2
pi_1 = 0.9; pi_2 = 0.1

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train2.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# Unbalanced Case 3 -------------------------------------------------------

# pi_1 > pi_2, pi_1 + pi_2 < 1.
# The propotion of each class in the Batch 1 are different. 
# The sample numbers of two batches are different.

## Case 3.1
pi_1 = 0.55; pi_2 = 0.35

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train3.1 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

## Case 3.2
pi_1 = 0.4; pi_2 = 0.1

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}


batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train3.2 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()


# None batch effect case --------------------------------------------------

# pi_1 = pi_2 = 1
# All samples come from the same batch

pi_1 = 1; pi_2 = 1

n_11 = 40*pi_1 # Samples in Class 1 from Batch 1
n_12 = 40 - n_11 # Samples in Class 1 from Batch 2
n_21 = 40*pi_2 # Samples in Class 2 from Batch 1
n_22 = 40 - n_21 # Samples in Class 2 from Batch 2

M = matrix(nrow = n, ncol = t4) # In the data matrix M, samples are in the rows and feasures are in the columns. 
colnames(M) = paste0('X', 1:t4)

if(t1 > 0){
  for(i in 1 : t1){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + rnorm(n, 0 ,1)
  }
}

if(t2 > t1){
  for(i in (t1 + 1) : t2){
    M[,i] = c(rep(2, 40), rep(-2, 40)) + 
      c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t3 > t2){
  for(i in (t2 + 1) : t3){
    M[,i] = c(rnorm(n_11, 2, 1), rnorm(n_12, -2, 1), rnorm(n_21, 2, 1), rnorm(n_22, -2, 1)) + 
      rnorm(n, 0 ,1)
  }
}

if(t4 > t3){
  for(i in (t3 + 1) : t4){
    M[,i] = rnorm(n, 0 ,1)
  }
}

batch = rep(c('Batch1', 'Batch2', 'Batch1', 'Batch2'), c(n_11, n_12, n_21, n_22) ) 
train4 = data.frame(y, M, batch, pi_1, pi_2)

## Visualize the training data 

Edata = t(M)
Mdata = melt(Edata)

ggplot(Mdata, aes(x = Var2, y = Var1, fill = value)) +
  labs(x = "Sample", y = "Gene", fill = "Value") +
  geom_tile() + 
  scale_fill_gradient2()

# save the data set -------------------------------------------------------

save(train0, train1.1, train1.2, train2.1, train2.2, 
     train3.1, train3.2, train4, fe, file = "rdata/mean_hetero_example.RData")
