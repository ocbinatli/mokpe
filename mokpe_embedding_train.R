mokpe_embedding_train <- function(K_c, K_x, K_z, parameters) {
  
  start_time <- Sys.time()
  N_x <- dim(K_x)[1]
  N_z <- dim(K_z)[1]
  R <- parameters$R
  sigma_e <- parameters$sigma_e
  lambda_c <- parameters$lambda_c
  lambda_x <- parameters$lambda_x
  lambda_z <- parameters$lambda_z
  iteration <- parameters$iteration
  tol <- parameters$tol
  seed <- parameters$seed
  set.seed(seed)
  
  indices_c <- which(is.na(K_c) == FALSE, arr.ind = T)
  count_c <- length(indices_c) / 2
  indices_x <- which(is.na(K_x) == FALSE, arr.ind = T)
  count_x <- length(indices_x) / 2
  indices_z <- which(is.na(K_z) == FALSE, arr.ind = T)
  count_z <- length(indices_z) / 2
  
  P_x <- rnorm(N_x * R)
  P_z <- rnorm(N_z * R)
  dim(P_x) <- c(N_x, R)
  dim(P_z) <- c(N_z, R)
  E_x <- project_to_stiefel_manifold(P_x)
  E_z <- project_to_stiefel_manifold(P_z)
  DE_c <- (as.matrix(pdist(E_x, E_z)))^2
  KE_c <- exp(-DE_c / sigma_e^2)
  
  DE_x <- (as.matrix(pdist(E_x, E_x)))^2
  KE_x <- exp(-DE_x / sigma_e^2)
  DE_z <- (as.matrix(pdist(E_z, E_z)))^2
  KE_z <- exp(-DE_z / sigma_e^2)
  objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
  objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
  objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
  objective <- lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  
  iter <- 0
  current_time <- Sys.time()
  # elapsing time in seconds
  timediff_sec <- as.double(difftime(current_time, start_time, units = "secs"))
  print(paste("iteration=", iter, "objective=", objective, "elapsed time=", timediff_sec))
  
  #---X
  tx <- function(x) {
    E_x = matrix(x[1:(dim(P_x)[1] * dim(P_x)[2])], dim(P_x)[1], dim(P_x)[2])
  }
  
  fx <- function(x) {
    E_x <- tx(x)
    
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    KE_c <- exp(-DE_c / sigma_e^2)
    DE_x <- (as.matrix(pdist(E_x, E_x)))^2
    KE_x <- exp(-DE_x / sigma_e^2)
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
    
    lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  }  
  
  gx <- function(x) {
    E_x <- tx(x)
    
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    KE_c <- exp(-DE_c / sigma_e^2)
    DE_x <- (as.matrix(pdist(E_x, E_x)))^2
    KE_x <- exp(-DE_x / sigma_e^2)
    
    E_x_gradient <- matrix(0L, nrow = N_x, ncol = R) 
    for (i in 1:N_x) {
      j <- !is.na(K_c[i,])
      if (sum(j) == 1) {
        E_x_gradient[i,] <- c(E_x_gradient[i,] - 4 * lambda_c * (as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_x[i,] - t(E_z[j,])) ) / sigma_e^2 / count_c)
      } else {
        E_x_gradient[i,] <- E_x_gradient[i,] - 4 * lambda_c * colSums(as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_x[i,] - t(E_z[j,])) ) / sigma_e^2 / count_c
      } 
    }
    
    if (lambda_x != 0) {
      for (i in 1:N_x) {
        j <- !is.na(K_x[i,])
        if (sum(j) == 1) {
          E_x_gradient[i,] <- c(E_x_gradient[i,] - 8 * lambda_x * (as.vector(t(KE_x[i, j]) * t(t(KE_x[i, j] - t(K_x[i, j]))) ) * t(E_x[i,] - t(E_x[j,])) ) / sigma_e^2 / count_x)
        } else {
          E_x_gradient[i,] <- E_x_gradient[i,] - 8 * lambda_x * colSums(as.vector(t(KE_x[i, j]) * t(t(KE_x[i, j] - t(K_x[i, j]))) ) * t(E_x[i,] - t(E_x[j,])) )  / sigma_e^2 / count_x
        }
      }
    }
    E_x_gradient
  }
  
  #---Z
  tz <- function(x) {
    E_z = matrix(x[1:(dim(P_z)[1] * dim(P_z)[2])], dim(P_z)[1], dim(P_z)[2])
  }
  
  fz <- function(x) {
    E_z <- tz(x)
    
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    KE_c <- exp(-DE_c / sigma_e^2)
    DE_z <- (as.matrix(pdist(E_z, E_z)))^2
    KE_z <- exp(-DE_z / sigma_e^2)
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
    
    lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  }
  
  gz <- function(x) {
    E_z <- tz(x)
    
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    KE_c <- exp(-DE_c / sigma_e^2)
    DE_z <- (as.matrix(pdist(E_z, E_z)))^2
    KE_z <- exp(-DE_z / sigma_e^2)
    
    E_z_gradient <- matrix(0L, nrow = N_z, ncol = R) 
    for (j in 1:N_z) {
      i <- !is.na(K_c[, j])
      if (sum(i) == 1) {
        E_z_gradient[j,] <- c(E_z_gradient[j,] - 4 * lambda_c * (as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_z[j,] - t(E_x[i,])) ) / sigma_e^2 / count_c)
      } else {
        E_z_gradient[j,] <- E_z_gradient[j,] - 4 * lambda_c * colSums(as.vector(t(KE_c[i, j]) * t(t(KE_c[i, j] - t(K_c[i, j]))) ) * t(E_z[j,] - t(E_x[i,])) ) / sigma_e^2 / count_c
      }
    }
    
    if (lambda_z != 0){
      for (i in 1:N_z){
        j <- !is.na(K_z[i,])
        if (sum(j) == 1) {
          E_z_gradient[i,] <- c(E_z_gradient[i,] - 8 * lambda_z * (as.vector(t(KE_z[i, j]) * t(t(KE_z[i, j] - t(K_z[i, j]))) ) * t(E_z[i,] - t(E_z[j,])) ) / sigma_e^2 / count_z)
        } else {
          E_z_gradient[i,] <- E_z_gradient[i,] - 8 * lambda_z * colSums(as.vector(t(KE_z[i, j]) * t(t(KE_z[i, j] - t(K_z[i, j]))) ) * t(E_z[i,] - t(E_z[j,])) ) / sigma_e^2 / count_z
        }
      }
    }
    E_z_gradient
  }
  
  #---eta
  teta <- function(x) {
    sigma_e = x[1:1]
  }
  
  feta <- function(x) {
    sigma_e <- teta(x)
    
    KE_c <- exp(-DE_c / exp(sigma_e)^2)
    KE_x <- exp(-DE_x / exp(sigma_e)^2)
    KE_z <- exp(-DE_z / exp(sigma_e)^2)
    
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
    objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
    
    lambda_c * objective_c + lambda_x * objective_x + lambda_z * objective_z
  }
  
  geta <- function(x) {
    sigma_e <- teta(x)
    
    KE_c <- exp(-DE_c / exp(sigma_e)^2)
    KE_x <- exp(-DE_x / exp(sigma_e)^2)
    KE_z <- exp(-DE_z / exp(sigma_e)^2)
    
    eta_gradient <- 0
    eta_gradient <- eta_gradient - 4 * lambda_c * sum(KE_c[indices_c] * (KE_c[indices_c] - K_c[indices_c]) * DE_c[indices_c]) / exp(sigma_e)^2 / count_c
    eta_gradient <- eta_gradient - 4 * lambda_x * sum(KE_x[indices_x] * (KE_x[indices_x] - K_x[indices_x]) * DE_x[indices_x]) / exp(sigma_e)^2 / count_x
    eta_gradient <- eta_gradient - 4 * lambda_z * sum(KE_z[indices_z] * (KE_z[indices_z] - K_z[indices_z]) * DE_z[indices_z]) / exp(sigma_e)^2 / count_z
    -1 * eta_gradient
  }
  
  mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  probx <- new(mod$RProblem, fx, gx)
  probz <- new(mod$RProblem, fz, gz)
  probeta <- new(mod$RProblem, feta, geta) 
  
  manix.defn <- get.stiefel.defn(dim(P_x)[1], dim(P_x)[2])
  maniz.defn <- get.stiefel.defn(dim(P_z)[1], dim(P_z)[2])
  manieta.defn <- get.euclidean.defn(1, 1)
  
  mani.params <- get.manifold.params()
  
  solver.params <- get.solver.params(DEBUG = 0)
  solver.params_eta <- get.solver.params(DEBUG = 0)
  
  # set the optimization algorithm
  method_algo <- "LRBFGS"	
  objective_c_all = c(objective_c)
  objective_x_all = c(objective_x)
  objective_z_all = c(objective_z)
  objective_all = c(objective)
  
  for (iter in 1:iteration) {
    #---X
    x0 <- as.numeric(E_x)
    resx <- manifold.optim(probx, mani.defn = manix.defn, method = method_algo,
                           mani.params = mani.params, solver.params = solver.params, x0 = x0)
    E_x <- tx(resx$xopt)
    DE_x <- (as.matrix(pdist(E_x, E_x)))^2
    KE_x <- exp(-DE_x / sigma_e^2)
    objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
    
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    KE_c <- exp(-DE_c / sigma_e^2)
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    objective_c_all <- c(objective_c_all, objective_c)
    objective_x_all <- c(objective_x_all, objective_x)
    objective_z_all <- c(objective_z_all, objective_z)
    objective_all <- c(objective_all, resx$fval)
    
    current_time <- Sys.time()
    timediff_sec_new <- as.double(difftime(current_time, start_time, units = "secs"))
    timediff_sec <- c(timediff_sec, timediff_sec_new)
    print(paste("iteration=", iter, "objective=", resx$fval, "elapsed time=", timediff_sec_new))
    
    #---Z
    z0 <- as.numeric(E_z)
    resz <- manifold.optim(probz, mani.defn = maniz.defn, method = method_algo,
                           mani.params = mani.params, solver.params = solver.params, x0 = z0)
    E_z <- tz(resz$xopt)
    DE_c <- (as.matrix(pdist(E_x, E_z)))^2
    DE_z <- (as.matrix(pdist(E_z, E_z)))^2
    
    KE_z <- exp(-DE_z / sigma_e^2)
    objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
    KE_c <- exp(-DE_c / sigma_e^2)
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    objective_c_all <- c(objective_c_all, objective_c)
    objective_x_all <- c(objective_x_all, objective_x)
    objective_z_all <- c(objective_z_all, objective_z)
    objective_all <- c(objective_all, resz$fval)
    
    current_time <- Sys.time()
    timediff_sec_new <- as.double(difftime(current_time, start_time, units = "secs"))
    timediff_sec <- c(timediff_sec, timediff_sec_new)
    print(paste("iteration=", iter, "objective=", resz$fval, "elapsed time=", timediff_sec_new))
    
    #---eta
    reseta <- manifold.optim(probeta, mani.defn = manieta.defn, method = method_algo,
                             mani.params = mani.params, solver.params = solver.params_eta, x0 = log(sigma_e))
    
    sigma_e <- exp(teta(reseta$xopt))
    KE_z <- exp(-DE_z / sigma_e^2)
    objective_z <- (sum((KE_z[indices_z] - K_z[indices_z])^2)) / count_z
    
    KE_c <- exp(-DE_c / sigma_e^2)
    objective_c <- (sum((KE_c[indices_c] - K_c[indices_c])^2)) / count_c
    KE_x <- exp(-DE_x / sigma_e^2)
    objective_x <- (sum((KE_x[indices_x] - K_x[indices_x])^2)) / count_x
    objective_c_all <- c(objective_c_all, objective_c)
    objective_x_all <- c(objective_x_all, objective_x)
    objective_z_all <- c(objective_z_all, objective_z)
    objective_all <- c(objective_all, reseta$fval)
    
    current_time <- Sys.time()
    timediff_sec_new <- as.double(difftime(current_time, start_time, units = "secs"))
    timediff_sec <- c(timediff_sec, timediff_sec_new)
    
    print(paste("iteration=", iter, "objective=", reseta$fval, "elapsed time=", timediff_sec_new))
    if ((length(objective_all) > 10) && (tail(objective_all)[3] - tail(objective_all)[6] < tol)) {
      break
    }
  }
  # set the state list
  state <- list(E_x = E_x, E_z = E_z, sigma_e = sigma_e, objective_c = objective_c_all, objective_x = objective_x_all, objective_z = objective_z_all, objective = objective_all, time = timediff_sec)
}

project_to_stiefel_manifold <- function(Q) {
  s <- svd(Q, nu = nrow(Q), nv = ncol(Q))
  Q <- s$u %*% diag(nrow = dim(Q)[1], ncol = dim(Q)[2]) %*% t(s$v)
}
