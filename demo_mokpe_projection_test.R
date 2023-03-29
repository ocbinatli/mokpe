library(pdist)
library(pracma)
library(pROC)

pdist <- pdist::pdist

source("mokpe_projection_test.R")

# N_x x N_x drug compound structure similarity matrix
cd <- sprintf("./%s_admat_dgc.txt", network)
# N_z x N_z protein sequence similarity matrix
wdd <- sprintf("./%s_simmat_dc.txt", network)
# N_x x N_z adjacency matrix containing drug-target interactions
wdp <- sprintf("./%s_simmat_dg.txt", network)

Y <- t(as.matrix(read.table(cd, header = TRUE)))
Y[Y == 0] <- -1
Y[Y == 1] <- 1

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

# set the seed to select the leftout drugs for each replication
set.seed(???)
fold_indices <- fold_indices[randperm(N_x)]

# load state
load(file = ???)

prediction <- rep(list(list()), fold)

for (f in 1:fold) {
  train_indices <- which(fold_indices != f, arr.ind = TRUE)
  test_indices <- which(fold_indices == f, arr.ind = TRUE)
  # perform prediction
  prediction[[f]] <-
    mokpe_projection_test(K_x_feature[train_indices, test_indices],
                         K_z_feature, state[[f]])
}
save("prediction", file = ???)

result_auc <- zeros(1, fold)

for (f in 1:fold) {
  test_indices <- which(fold_indices == f, arr.ind = TRUE)
  target <- Y[test_indices,]
  pd <- (as.matrix(pdist(prediction[[f]]$E_x, prediction[[f]]$E_z))) ^ 2
  predicted <- exp(-pd / state[[f]]$sigma_e ^ 2)
  calculated_auc <-
    auc(roc(
      response = as.vector(target),
      predictor = as.vector(predicted),
      levels = c(1,-1)
    ))
  if (is.na(calculated_auc) == FALSE) {
    result_auc[f] <- calculated_auc
  } else {
    result_auc[f] <- 1
  }
}

save("result_auc", file = ???)