#' Test of local structure for spatio-temporal point processes on linear networks
#'
#' @param X Background spatio-temporal point pattern on a linear network. Must be an stlnpp object.
#' @param Z Other spatio-temporal point pattern on a linear network. Must be an stlnpp object.
#' @param method Character string indicating which version of LISTA function to use: either "K" or "g".
#' If "K", the Local spatio-temporal K-function on a linear network is used to run the test.
#' If "g", the Local spatio-temporal pair correlation function is used.
#' @param k Number of permutations
#' @param parallel FALSE is default. If TRUE, parallelization is performed.
#' @param ncores Number of parallel threads. Optional.
#'
#' @import parallel
#'
#' @return A vector of p-values, one for each of the points in X.
#' @export
#'
#' @examples
#'
#'
#'
#'\dontrun{
#'#simulate two spatio-temporal point patterns on a linear network with different local structure
#'L0 = domain(chicago)
#'set.seed(4)
#'X0 = sim_ETASnet(cat = NULL, params = c(0.078915/2, 0.003696,  0.013362,  1.2,0.424466,  1.164793),
#'                          betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
#'                          xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300,
#'                          iprint = TRUE, covdiag = FALSE, covsim = FALSE, L = L0)$cat.sim
#'X = as.stlpp(X0$long, X0$lat, X0$time, L = L0)
#'lambda0 <- 20 / (volume(domain(chicago)) * (200-25))
#'set.seed(4)
#'Z = rpoistlpp(lambda = lambda0, a = 25, b = 200, L = L0)
#'
#'#compute the local test
#'set.seed(1); system.time(p <- test_local(X, Z, method = "K", k = 19, parallel = FALSE))
#'
#'#display the significant and non-significan points
#'plot(X[which(p <= 0.05)])
#'plot(X[-which(p <= 0.05)])
#'}
#'
#'
test_local <- function(X, Z, method = "K", k, parallel = FALSE, ncores = 4) {

  if (!(inherits(X,"stlpp")|inherits(Z,"stlpp")))
    stop("X and Z should be both from class stlpp")

  if(parallel) cl <- makeCluster(getOption("cl.cores", ncores))

  n_pointsX <- npoints(X)
  n_pointsZ <- npoints(Z)
  perm <- c(rep(1, n_pointsZ), rep(2, n_pointsX - 1))
  # Tiq <- matrix(NA, k, n_pointsX)
  p <- double(n_pointsX)
  L0 <- domain(X)
  res_test <- double(n_pointsX)
  res_test_d <- matrix(NA, k, n_pointsX)

  if (method == "K"){
    X_lista = STLKinhom_local_a(X, lambda = rep(npoints(X) / (volume(X$domain) * (X$time[2] - X$time[1])), npoints(X)))
  } else if (method == "g"){
    X_lista = STLginhom_local_a(X, lambda = rep(npoints(X) / (volume(X$domain) * (X$time[2] - X$time[1])), npoints(X)))
  } else {stop(" 'method' argument must be either \"g\" or \"K\ ")}

  # AVG <- rowMeans(X_lista, dims = 2)
  lista_0 <- array(NA, dim = c(dim(X_lista)[1:3], k))

  for (i in 1:n_pointsX){#punti del processo
    print(i)

    Q_n <- data.frame(rbind(cbind(Z$data, 1), cbind(X[-i]$data, 2)))

    if (parallel){
      clusterExport(cl = cl, c('permut_test', 'as.stlpp', 'STLKinhom_local_a', 'volume', 'as.lpp.stlpp', 'domain', 'npoints', 'pairdist.lpp',
                               'default.linnet.tolerance', 'countends'))
      lista_0 <- parSapply(cl, 1:k, function(q) permut_test(perm = perm, Q_n = Q_n, Xi = X[i], L0 = L0))
    }
    else {
      pb <- txtProgressBar(min = 0, max = k, style = 3)
      lista_0 <- sapply(1:k, function(q) {
        setTxtProgressBar(pb, q)
        permut_test(perm = perm, Q_n = Q_n, Xi = X[i], L0 = domain(X))
      }, simplify = "array")
      close(pb)
    }
    # AVG0 <- rowMeans(Q_lista, dims = 2)
    # Tiq <- sum((Q_lista[, , 1] - AVG0)^2)
    # Tiq
    # Ti <- sum((X_lista[, , i] - AVG) ^ 2)
    # p[i] <- mean(Tiq >= Ti)
  }

  mlistas <- apply(lista_0, c(1, 2, 3), mean, na.rm = TRUE) # Mean of the curves
  # varlistas <- apply(lista_0,c(1,2,3),var,na.rm=TRUE) # standar desviation

  for (i in 1:n_pointsX){
    lista_a <- X_lista[, , i]  #lista for the ith point
    lista_b <- mlistas[, , i]
    # lista_c<-varlistas[,,i]
    t2 <- sum((lista_a - lista_b) ^ 2)
    res_test[i] <- t2
    # t4<-sum((lista_a-lista_b)^2/lista_c,na.rm=T)
    for (j in 1:k){
      lista_a <- lista_0[, , i, j]
      t2_0 <- sum((lista_a - lista_b) ^ 2)
      # t4_0<-sum((lista_a-lista_b)^2/lista_c,na.rm=T)
      res_test_d[j,i] <- t2_0
    }
    p[i] <- mean(res_test_d[,i] >= res_test[i])
  }

  if (parallel) stopCluster(cl)

  return(p)

}
