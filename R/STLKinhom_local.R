#' Local inhomogeneous Spatio-temporal K functions on a linear network
#'
#' This function calculates the local inhomogeneous K-function for a spatio-temporal point processes on a linear network.
#' @param X a realisation of a spatio-temporal point processes on a linear networks.
#' @param lambda values of estimated intensity.
#' @param normalize normalization factor to be considered.
#' @param r values of argument r where K-function will be evaluated. optional.
#' @param t values of argument t where K-function will be evaluated. optional.
#' @param nxy pixel array dimensions. optional.
#' @keywords
#' @import spatstat stlnpp
#' @return A list of objects of class sumstlpp.
#' @export
#' @examples
#'
#'\dontrun{
#'set.seed(10)
#'X <- rpoistlpp(.2, a = 0, b = 5, L = easynet)
#'X
#'lambda <- density(X, at = "points")
#'k <- STLKinhom_local(X, lambda = lambda, normalize = TRUE)
#'#select an individual point
#'j = 1
#'k[[j]]
#'#plot the lista function and compare it with its theoretical value
#'inhom <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Kinhom)
#'theo <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Ktheo)
#'diff <- list(x = k[[j]]$r, y = k[[j]]$t, z = k[[j]]$Kinhom - k[[j]]$Ktheo)
#'par(mfrow=c(1,3))
#'fields::image.plot(inhom, main= "Kinhom", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'fields::image.plot(theo, main = "Ktheo", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'fields::image.plot(diff, main = "Kinhom - Ktheo", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'}
#'
#'

STLKinhom_local = function(X, lambda = lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
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
  edgetl <- mtedge * ml * lamden

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
        kout <- out /edgetl
        kout <- kout[, k]#
        K[i, j] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
      }
    }
    K_local[[k]] <- K#
  }#

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1) * (tleng * trange) / (sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)#
    K_local <- lapply(K_local, "*" , appx)#

  }
  else {
    K_local <- lapply(K_local, FUN= function(K) (n - 1) * K / (trange * tleng) )#
  }
  K_theo <- list()#
  for(k in 1:n){#
    K_theo[[k]] <- matrix(expand.grid(r, t)[, 1] * expand.grid(r, t)[,2],ncol = nxy)#
  }#

  Kout <- list()#
  for(k in 1:n){#
    Kout[[k]] <- list(Kinhom = K_local[[k]], Ktheo = K_theo[[k]], r = r, t = t)#
  }#

  class(Kout[[k]]) <- c("sumstlpp")#
  attr(Kout[[k]], "nxy") <- nxy#
  return(Kout)
}


