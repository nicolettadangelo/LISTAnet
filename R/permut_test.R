#' Permutation test for a spatio- temporal point process on a linear network
#'
#' @param perm Number of permutations to use
#' @param Q_n Data frame containing the spatio-temporal locations of the points to be permuted
#' @param Xi The point of the original point pattern
#' @param L0 Linear network of the spatio-temporal point process to be permuted
#'
#' @return A vector of squared discrepancy values from the lista functions computed on the points of the pattern and their average
#'
#' @examples
permut_test <- function(perm, Q_n, Xi, L0){
  n_perm <- length(perm)
  n_points <- sum(perm == 2) + 1
  A <- sample(perm, size = n_perm, replace = FALSE, prob = NULL)
  cod <- which(A == 2)
  Q_perm <- Q_n[cod, 1:3]
  Q <- as.stlpp(c(Xi$data$x[[1]], Q_perm$x),
                c(Xi$data$y[[1]], Q_perm$y),
                c(Xi$data$t[[1]], Q_perm$t),
                L = L0)

  Q_dens <- rep(n_points / (volume(Q$domain) * (Q$time[2] - Q$time[1])), n_points)
  lista_0 <- STLKinhom_local_a(Q, lambda = Q_dens)
  lista_0
  # AVG0 <- rowMeans(Q_lista, dims = 2)
  # Tiq <- sum((Q_lista[, , 1] - AVG0)^2)
  # Tiq
}
