#' @importFrom MASS mvrnorm
#' @export
getQT1 <- function(n,p,theta,B)
{
  p_tilde <- min(n,p)
  output <- matrix(0, nrow = p_tilde, ncol = B)

  for (b in 1:B)
  {
    Sigma <- diag(rgamma(n = p, shape = theta, rate = theta))
    X <- mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
    if (n < p) S <- 1 / n * tcrossprod(X) ## tcrossprod(X) = XX'
    else S <- 1 / n * crossprod(X) ## crossprod(X) = X'X
    eigen.S <- eigen(S)
    output[,b] <- eigen.S$values
  }
  rowMeans(output)

}
