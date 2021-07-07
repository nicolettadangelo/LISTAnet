#' Simulation of a spatio-temporal LGCP (Log-Gaussian Cox Process) model on a linear network
#'
#' This function simulates a spatio-temporal LGCP (Log-Gaussian Cox Process) model on a linear network.
#' @param lambda
#' @param a
#' @param b
#' @param L
#' @param check
#' @param lmax
#' @param nsim
#' @param sigma
#' @param phi
#' @param theta
#' @param nmark
#' @keywords
#' @import spatstat stlnpp
#' @return an object of class stlpp if nsim=1, otherwise a list of objects of class stlpp
#' @export
#' @examples
#'
#'
#'\dontrun{
#'L0 = domain(chicago)
#'aa = 64 / ((1 - exp( - 4)) * (1 - exp( - 2)))
#'lambda0 = function(x, y, t) {exp(log(aa) - 4 * y - 2 * t)}
#'X = sim_LGCPnet(lambda0, a = 0, b = 4, nmark = 1, sigma = 2, phi = 0.5, theta = 1,L = L0, lmax = 10)
#'X
#'plot(X)
#'}
#'
#'




sim_LGCPnet = function (lambda = lambda0,
                 a = 0,
                 b = 3,
                 L = L0,
                 check = FALSE, lmax = NULL, nsim = 1,
                 sigma = 2,
                 phi = 0.5,
                 theta = 2,
                 nmark = 1)
{

  if (!inherits(L, "linnet"))
    stop("L should be a linear network")

  if (a >= b)
    stop("lower bound must be smaller than upper bound")

  if (!is.numeric(lambda) & !is.function(lambda))
    stop(" lambda should be a number or a function")

  if (nsim > 1) {
    out <- list()
    for (i in 1:nsim) {
      out[[i]] <- rpoislpp(lambda, a = a, b = b, L = L,
                           nsim = 1)
    }
    return(out)
  }
  if (is.numeric(lambda)) {
    n <- rpois(1, lambda * volume(L) * (b - a))
    X <- runiflpp(n, L)
    t <- runif(npoints(X), a, b)
    stlpp <- data.frame(x = X$data$x, y = X$data$y, t)
  }
  else {
    if (is.null(lmax)) {
      Llines <- as.psp(L)
      linemask <- as.mask.psp(Llines,dimyx=4)
      lineimage <- as.im(linemask)
      xx <- raster.x(linemask)
      yy <- raster.y(linemask)
      mm <- linemask$m
      xx <- as.vector(xx[mm])
      yy <- as.vector(yy[mm])
      pixelcentres <- ppp(xx, yy, window = as.rectangle(linemask),
                          check = check)
      pixelcentres <- unique.ppp(pixelcentres)
      pixdf <- data.frame(xc = xx, yc = yy)
      p2s <- project2segment(pixelcentres, Llines)
      projloc <- as.data.frame(p2s$Xproj)
      projmap <- as.data.frame(p2s[c("mapXY", "tp")])
      projdata <- cbind(pixdf, projloc, projmap)
      gridx <- p2s$Xproj$x
      gridy <- p2s$Xproj$y
      df <- data.frame(gridx, gridy)
      df <- df[!duplicated(df), ]
      grid <- lpp(df, L)
      t0 <- purrr::rdunif(npoints(grid), a, b)
      prov=cbind(grid$data$x,grid$data$y,t0)
      colnames(prov)=c("x","y","t")
      prov=as.data.frame(prov)
      #grf
      t <- a:b
      sim <- CompRandFld::RFsim(grid$data$x, grid$data$y, t, corrmodel = "exp_exp", grid = FALSE,
                param = list(nugget = 0,mean = - sigma / 2, scale_s = phi, scale_t = theta, sill = 1))#sim$data
      newdata <- t(sim$data)
      newdata2 <- ks::vec(newdata)
      newdata3 <- cbind(newdata2, rep(sim$coordx, each = dim(sim$data)[1]), rep(sim$coordy, each = dim(sim$data)[1]), rep(sim$coordt, dim(sim$data)[2]))
      colnames(newdata3) <- c("int", "x0", "y0", "t0")
      newdata3 <- as.data.frame(newdata3)
      newdata3$id <- paste0(newdata3$x0, newdata3$y0, newdata3$t0)
      prov$id <- paste0(prov$x, prov$y, prov$t)
      newdata4 <- full_join(newdata3, prov)
      newdata5 <- subset(newdata4, !is.na(t))
      grf0 <- newdata5$int
      #int
      mark.prova <- sample(as.factor(1:nmark), npoints(grid), replace = TRUE)
      int0 <- lambda(grid$data$x, grid$data$y, t0)
      int <- int0 * exp(grf0)
      int <- na.omit(int)
      lmax<- max(int)
    }
    mean <- lmax * volume(L) * (b - a)
    n <- rpois(1, mean)
    unipoint <- runiflpp(n, L)
    t0 <- purrr::rdunif(n, a, b)
    sim <- CompRandFld::RFsim(unipoint$data$x, unipoint$data$y, t, corrmodel = "exp_exp", grid = FALSE,
              param = list(nugget = 0, mean = - sigma / 2, scale_s = phi, scale_t = theta, sill = 1))#sim$data
    hlpp <- cbind(unipoint$data$x, unipoint$data$y, t0)
    colnames(hlpp) <- c("x", "y", "t")
    hlpp <- as.data.frame(hlpp)
    newdata <- t(sim$data)
    newdata2 <- ks::vec(newdata)
    newdata3 <- cbind(newdata2, rep(sim$coordx, each = dim(sim$data)[1]), rep(sim$coordy, each = dim(sim$data)[1]), rep(sim$coordt, dim(sim$data)[2]))
    colnames(newdata3) <- c("int", "x0", "y0", "t0")
    newdata3 <- as.data.frame(newdata3)
    newdata3$id <- paste0(newdata3$x0, newdata3$y0, newdata3$t0)
    hlpp$id <- paste0(hlpp$x, hlpp$y, hlpp$t)
    newdata4 <- full_join(newdata3, hlpp)
    newdata5 <- subset(newdata4, !is.na(t))
    grf0 <- newdata5$int
    #int
    mark.prova2 <- sample(as.factor(1:nmark), length(hlpp[, 2]), replace = TRUE)
    int0 <- lambda(unipoint$data$x, unipoint$data$y, t0)
    int <- int0 * exp(grf0)
    int <- na.omit(int)
    prob <- int / max(int)
    #prob=int/lmax
    if (check) {
      if (any(prob < 0))
        warning("Negative values of lambda obtained")
      if (any(prob > 1))
        warning("lmax is not an upper bound for lambda")
    }
    u <- runif(length(hlpp[, 1]))
    retain <- (u <= prob)
    stlpp <- hlpp[retain, ]
    stlpp <- data.frame(x = stlpp[, 1], y = stlpp[, 2], t = stlpp[,
                                                                  3])
  }
  out <- ppx(data = stlpp, domain = L, coord.type = c("s",
                                                      "s", "t"))
  print(npoints(out))
  class(out) <- c("stlpp", "ppx")
  out$time <- c(a, b)
  out$lmax <- lmax
  return(out)
}

