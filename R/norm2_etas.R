
#'
#'
norm2_etas = function(n = 1, d = 1, q = 2, x0, y0){
  theta   = runif(n, max = 2 * pi)
  U       = runif(n)
  R       = sqrt(d * (U ^ (1 / (1 - q)) - 1))
  return(cbind(x0 + R * cos(theta), y0 + R * sin(theta)))
}
