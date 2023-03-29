library(ManifoldOptim)
library(pdist)
pdist <- pdist::pdist

source("mokpe_embedding_train.R")

# N_x x N_x drug compound structure similarity matrix
wdd <- sprintf("./%s_simmat_dc.txt", network)
# N_z x N_z protein sequence similarity matrix
wdp <- sprintf("./%s_simmat_dg.txt", network)
# N_x x N_z adjacency matrix containing drug-target interactions
cd <- sprintf("./%s_admat_dgc.txt", network)

# set the seed to initalize random variables
seed <- ???

K_c <- t(as.matrix(read.table(cd, header = TRUE)))
K_c[K_c == 0] <- NA
K_c <- K_c * 0.9

K_x <- t(as.matrix(read.table(wdd, header = TRUE)))
K_x <- (K_x + t(K_x)) / 2

K_z <- as.matrix(read.table(wdp, header = TRUE))
K_z <- (K_z + t(K_z)) / 2

# set the maximum number for (outer) iterations
iteration <- 10000

# set the stopping criteria for (outer) iterations
tol <- 1e-6

# set the regularization parameter for drug-target interactions
lambda_c <- 1.0

# set the regularization parameters for drug-drug, target-target similarities
lambda_x <- 0.1
lambda_z <- 0.1

# set the subspace dimensionality
R <- 2

# set the kernel width used in the subspace
sigma_e <- sqrt(R)

parameters <- data.frame(iteration, tol, lambda_c, lambda_x, lambda_z, R, seed, sigma_e)

# perform training
state <- mokpe_embedding_train(K_c, K_x, K_z, parameters)

# display the embeddings
# print(state$E_x)
# print(state$E_z)

save("state", file = ???)