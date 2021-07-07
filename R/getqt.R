#' @export
getQT1 <- function(n,p,theta,B)
{
  p_tilde <- min(n,p)
  output <- matrix(0, nrow = p_tilde, ncol = B)

  for (b in 1:B)
  {
    ## this is Sigma^1/2
    ## this approach should be faster than using mvrnorm
    Sigma_half <- diag(sqrt(rgamma(n = p, shape = theta, rate = theta)))
    X_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
    X <- X_raw %*% Sigma_half
    # if (n < p) S <- 1 / n * tcrossprod(X) ## tcrossprod(X) = XX'
    # else S <- 1 / n * crossprod(X) ## crossprod(X) = X'X
    # eigen.S <- eigen(S)
    # output[,b] <- eigen.S$values
    output[,b] <- 1/n * svd(X, nu = 0, nv = 0)$d^2 ## this avoids computing the eigenvectors
  }
  rowMeans(output)

}
