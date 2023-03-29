mokpe_projection_test <- function(X, Z, state) {
  E_x <- t(X) %*% state$Q_x
  E_z <- t(Z) %*% state$Q_z
  prediction <- list(E_x = E_x, E_z = E_z)
}