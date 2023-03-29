library(pdist)
library(pracma)
library(ManifoldOptim)

pdist <- pdist::pdist

source("mokpe_projection_train.R")

# N_x x N_x drug compound structure similarity matrix
wdd <- sprintf("./%s_simmat_dc.txt", network)
# N_z x N_z protein sequence similarity matrix
wdp <- sprintf("./%s_simmat_dg.txt", network)
# set the seed to determine the drugs that will be leftout for each replication
set.seed(???)

K_x <- as.matrix(read.table(wdd, header = TRUE))
K_x <- (K_x + t(K_x)) / 2

K_z <- as.matrix(read.table(wdp, header = TRUE))
K_z <- (K_z + t(K_z)) / 2

K_x_feature <- K_x
K_z_feature <- K_z

N_x <- dim(K_x)[1]
for (i in 1:N_x) {
  K_x[i, i] <- NA
}

N_z <- dim(K_z)[1]
for (i in 1:N_z) {
  K_z[i, i] <- NA
}

fold <- 10
fold_indices <- 1:N_x
fold_indices <- mod(fold_indices + fold - 1, fold) + 1
fold_indices <- fold_indices[randperm(N_x)]

state <- rep(list(list()), fold)

for(f in 1:fold) {
  # N_x x N_z adjacency matrix containing drug-target interactions
  cd <- sprintf("./%s_admat_dgc.txt", network)
  
  Y <- t(as.matrix(read.table(cd, header = TRUE)))
  Y[Y == 0] <- NA
  Y <- Y * 0.9  
  
  # set the maximum number of (outer) iterations
  iteration <- 20
  
  # set the regularization parameter for drug-target interactions
  lambda_c <- 1.0
  
  # set the regularization parameters for drug-drug, target-target similarities
  lambda_x <- 0.1
  lambda_z <- 0.1
  
  # set the seed to initalize random variables for each fold
  seed <- ???
    
  # set the kernel width used in the subspace
  sigma_e <- sqrt(R)
  
  # set the parameters data frame
  parameters <- data.frame(iteration, lambda_c, lambda_x, lambda_z, R, seed, sigma_e)
  
  # perform training
  state[[f]] <- mokpe_projection_train(Y[fold_indices != f,], K_x[fold_indices != f, fold_indices != f], 
                                       K_z, K_x_feature[fold_indices != f, fold_indices != f], K_z_feature, parameters)
}
save("state", file = ???)
