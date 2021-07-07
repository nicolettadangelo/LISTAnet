
#' Title
#'
#' @param X
#' @param lambda
#' @param normalize
#' @param r
#' @param t
#' @param nxy
#'
#' @return
#'
#' @examples
STLKinhom_local_a = function(X, lambda = lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")

  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))
  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  diag(lamden) <- 1
  edgetl <- mtedge* ml * lamden

  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy), maxs, by = (maxs - (maxs/nxy))/(nxy - 1))
  if (is.null(t))
    t <- seq((maxt/nxy), maxt, by = (maxt - (maxt/nxy))/(nxy - 1))

  K_local <- list()#
  for(k in 1:n){#
    K <- matrix(NA, nrow = nxy, ncol = nxy)
    for (i in 1:length(r)) {
      for (j in 1:length(t)) {
        out <- (sdist <= r[i]) * (tdist <= t[j])
        diag(out) <- 0
        kout <- out / edgetl
        kout <- kout[, k]
        K[i, j] <- sum(kout[is.finite(kout)])
      }
    }
    K_local[[k]] <- K#
  }#

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1) * (tleng * trange) / (sum(revrho[lower.tri(revrho, diag = FALSE)]) * 2)#
    K_local <- lapply(K_local, "*" , appx)#

  }
  else {
    K_local <- lapply(K_local, FUN= function(K) (n - 1) * K / (trange * tleng) )#
  }

  arr <- array(unlist(K_local), c(10, 10, n))

  return(arr)
}
