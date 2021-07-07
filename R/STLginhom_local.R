#' Local inhomogeneous Spatio-temporal pair correlation functions on a linear network
#'
#' This function calculates the local inhomogeneous pair correlation function for a spatio-temporal point processes on a linear network.
#' @param X a realisation of a spatio-temporal point processes on a linear networks.
#' @param lambda values of estimated intensity.
#' @param normalize normalization factor to be considered.
#' @param r values of argument r where pair correlation function will be evaluated. optional.
#' @param t values of argument t where pair correlation function will be evaluated. optional.
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
#'g <- STLginhom_local(X, lambda = lambda, normalize = TRUE)
#'#select an individual point
#'j = 1
#'g[[j]]
#'#plot the lista function and compare it with its theoretical value
#'inhom <- list(x = g[[j]]$r, y = g[[j]]$t, z = g[[j]]$ginhom)
#'theo <- list(x = g[[j]]$r, y = g[[j]]$t, z = g[[j]]$gtheo)
#'diff <- list(x = g[[j]]$r, y = g[[j]]$t, z = g[[j]]$ginhom - g[[j]]$gtheo)
#'par(mfrow=c(1,3))
#'fields::image.plot(inhom, main = "ginhom", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'fields::image.plot(theo, main = "gtheo", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'fields::image.plot(diff, main = "ginhom - gtheo", col = hcl.colors(12, "YlOrRd", rev = FALSE), xlab = "Spatial distance", ylab = "Temporal distance")
#'}
#'
#'
#'
#'

STLginhom_local = function(X, lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(X, "stlpp"))
    stop("X should be from class stlpp")
  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(Y)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- pairdist(Y)
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
  edgetl <- mtedge * ml * lamden
  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if (is.null(t))
    t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))

  g_local = list()#
  for(k in 1:n){#

    g <- matrix(NA, nrow = nxy, ncol = nxy)
    no <- sdist == 0&tdist == 0|sdist == Inf|sdist > maxs|tdist > maxt
    bwl <- bw.nrd0(as.numeric(sdist[!no]))#bandwidth
    bwt <- bw.nrd0(as.numeric(tdist[!no]))#bandwidth

    for (i in 1:length(r)) {
      for (j in 1:length(t)) {

        outl <- dkernel(as.numeric(sdist[!no] - r[i]), sd = bwl)
        outt <- dkernel(as.numeric(tdist[!no] - t[j]), sd = bwt)
        g1 <- outl * outt / (edgetl[!no])
        no2 <- no#
        no2[no] <- NA#
        no2[!no] <- g1#
        no2 <- no2[, k]#chech this
        g[i, j] <- sum(no2[!is.na(no2) & !is.infinite(no2)])

      }
    }
    g_local[[k]] <- g#
  }#

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1)*(tleng * trange)/(sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)
    gval <- lapply(g_local, "*" ,appx)#
  }
  else {
    gval <- lapply(g_local, FUN= function(K) (n-1) * K / (trange * tleng) )#

  }

  g_theo <- list()#
  for(k in 1:n){#
    g_theo[[k]] <- matrix(rep(1, length(t) * length(r)), ncol = nxy)#
  }#

  gout <- list()#
  for(k in 1:n){
    gout[[k]] <- list(ginhom = gval[[k]], gtheo = g_theo[[k]], r = r, t = t)
  }
  class(gout[[k]]) <- c("sumstlpp")
  attr(gout[[k]], "nxy") <- nxy
  return(gout)
}

