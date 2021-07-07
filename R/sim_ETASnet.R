#' Simulation of a spatio-temporal ETAS (Epidemic Type Aftershock Sequence) model on a linear network
#'
#' This function simulates a spatio-temporal ETAS (Epidemic Type Aftershock Sequence) model on a linear network.
#' @param params A vector of parameters
#' @param betacov
#' @param m0
#' @param b
#' @param tmin Minimum value of time.
#' @param t.lag
#' @param xmin
#' @param xmax
#' @param ymin
#' @param ymax
#' @param iprint Default TRUE
#' @param covdiag Default FALSE
#' @param covsim Default FALSE
#' @param L linear network
#' @keywords
#' @import spatstat stlnpp
#' @return An ETAS catalogue on a linear network. It can be easily converted into an stlnpp object
#' @export
#' @examples
#'
#'\dontrun{
#'L0 = domain(chicago)
#'set.seed(95)
#'X0 = sim_ETASnet(cat = NULL, params = c(0.1293688525, 0.003696, 0.013362, 1.2,0.424466,  1.164793),
#'                          betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
#'                          xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300, iprint = TRUE,
#'                          covdiag = FALSE, covsim = FALSE, L = L0)$cat.sim
#'X <- as.stlpp(x = X0$long, y = X0$lat, t = X0$time,  L = L0)
#'X
#'plot(X)
#'}
#'
#'

sim_ETASnet <- function(cat=NULL
                                  , params = c(0.078915, 0.003696,  0.013362,  1.2,0.424466,  1.164793)
                                  , betacov = 0.39
                                  , m0 = 2.5
                                  , b = 1.0789
                                  , tmin = 0
                                  , t.lag = 200
                                  , xmin = 600
                                  , xmax = 2200
                                  , ymin = 4000
                                  , ymax = 5300
                                  , iprint = TRUE
                                  , covdiag = FALSE
                                  , covsim = FALSE
                                  , L){
  # require(spatstat)
  # require(stlnpp)
    ## version of 2019 Jajn 7th; covdiag usued

    xmin = L$window$xrange[1]#
    xmax = L$window$xrange[2]#
    ymin = L$window$yrange[1]#
    ymax = L$window$yrange[2]#

    mu	= params[1]
    k0  	= params[2]
    c   	= params[3]
    p   	= params[4]
    gamma   =0              ## obliged in this first version
    #        a   	= params[5] # to change because of covariates
    d       = params[5]
    q       = params[6]
    ncov    = length(betacov)
    #        covdiag =(ncov>1)
    covdiag =! covsim
    beta    = log(10) * b
    mm      = (ymax - ymin) / (xmax - xmin)
    qq      = ymin - mm * xmin
    #    print(n)
    cat.pois = NULL
    cat.sim = NULL
    cat.new = NULL
    tmax    = tmin + t.lag
    ak      = k0 * c ^ (1-p) / (p-1)
    sk      = (pi * d ^ (1 - q)) / (q - 1) ## added 10-10-2018; check zhuang 2015 (corssa)
    muback  = mu * (tmax - tmin)
    n0      = rpois(1, muback)
    nstart  = n0
    #print(c(n,n0))
    #print(c("nstart",nstart))
    if (nstart > 0) {

        lgen       = array(0, nstart)
        ind     = array(0, nstart)
        father  = array(0, nstart)

        ## simulation can start
        ## build the poisson catalog
        tback   = runif(n0, tmin, tmax)
        # xback   =runif(n0,xmin,xmax)#
        # yback   =runif(n0,ymin,ymax)#
        xback = runiflpp(n0, L)$data$x#
        yback = runiflpp(n0, L)$data$y#

        mback   = m0 + rexp(n0, rate = beta)
        zback   = array(0, n0)
        if (ncov == 1) {
            cov2back = array(0, n0)
        }
        else{            ## covdiag segment for cov2>1
            if (covsim)
            {
                cov2back = rnorm(n0) ^ 2
            }
            else
            {
                cov2back = abs(yback - mm * xback - qq) / (sqrt(1 + mm * mm))
            }
        }


        cat.new = cbind(tback, xback, yback, mback, zback, cov2back)
        cat.pois = cat.new

        cat.new             = cbind(cat.new, lgen, ind, father)
        colnames(cat.new)   = c("time", "long", "lat", "magn1", "z", "cov2", "lgen", "ind", "father")
        cat.new             = as.data.frame(cat.new)
        cat.new = cat.new[order(cat.new$time), ]
        i = 0
        while(!prod(cat.new$ind)){
            ## try to use i-th event as a root for aftershoks
            i   = i + 1
            if (iprint & ((i %% 100) == 0)) print(c(i, cat.new$time[i]))
            cat.new$ind[i]  = TRUE
            pred            = (cat.new$magn1[i] - m0) * betacov[1]
            if (ncov > 1) pred = pred + cat.new$cov2[i] * betacov[2]
            nexpected       = ak * sk * exp(pred)
            ni              = rpois(1, nexpected)
            if (ni > 0) {
                xy              = norm2_etas(ni, d, q, cat.new$long[i], cat.new$lat[i])
                ## instead of t=rpareto(ni,lambda=p-1,a=c)-c+cat.new$time[i], inversion of F(t) has been used
                t             = c * runif(ni) ^ ( - 1 / (p - 1)) - c + cat.new$time[i]
                #            print("---------------------------------------------------------")
                ind.txy = (t > tmin)&(t < tmax)&(xy[,1] > xmin)&(xy[,1] < xmax)&(xy[,2] > ymin)&(xy[,2] < ymax)
                #            print(cat.new[i,1:4])
                #            print(cbind(t,xy,ind.txy))
                ntxy = sum(ind.txy)
                if (ntxy > 0) {
                    m1      = m0 + rexp(ntxy, rate = beta)
                    z1      = array(0, ntxy)
                    cov2    = array(0, ntxy)
                    lgen    = array(1, ntxy) + cat.new$lgen[i]
                    ind     = array(0, ntxy)
                    father  = array(i, ntxy)
                    if (ncov > 1)
                    {            ## interpolation segment for cov2

                        if (covsim)
                        {
                            cov2 = rnorm(ntxy) ^ 2
                        }
                        else
                        {
                            xo      = xy[ind.txy,1]
                            yo      = xy[ind.txy,2]
                            cov2    = abs(yo - mm * xo - qq) / (sqrt(1 + mm * mm))
                        }


                    }
                    cat.son = cbind(t[ind.txy], xy[ind.txy, 1], xy[ind.txy, 2], m1, z1, cov2, lgen, ind, father)
                    colnames(cat.son)   = c("time", "long", "lat", "magn1", "z", "cov2", "lgen", "ind", "father")
                    cat.son = as.data.frame(cat.son)
                    #               print(cat.son)
                    cat.new = rbind(cat.new, cat.son)
                    cat.new = cat.new[order(cat.new$time), ]

                }
            }
        }
    }
    nson = nrow(cat.new) - n0
    return(list(cat=cat, cat.pois = cat.pois, cat.sim = cat.new, n0 = n0, nson = nson))
}


