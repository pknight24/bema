loss <- function(theta, sample_eigenvalues, n, p, alpha, B)
{

  simulated_quantiles <- getQT1(n = n, p = p, theta = theta, B = B)

  bulk_indices <- floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))

  sigma2_hat <- lm(sample_eigenvalues[bulk_indices] ~ simulated_quantiles[bulk_indices] - 1)$coef[[1]]
  error <- sum( (sample_eigenvalues[bulk_indices] - sigma2_hat * simulated_quantiles[bulk_indices])^2 )
  return(error)
}

#' @importFrom stats optimize
#' @export
step1 <- function(sample_eigenvalues, n, p, alpha, theta_interval = c(0.1, 50), B = 20)
{
  p_tilde <- min(n,p)
  l_hat <- sample_eigenvalues
  opt_out <- optimize(loss, interval = theta_interval, sample_eigenvalues = sample_eigenvalues,
                      n = n, p = p, alpha = alpha, B = B)
  theta_hat <- opt_out$minimum
  simulated_quantiles <- getQT1(n = n, p = p, theta = theta_hat, B = B)
  bulk_indices <- floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  sigma2_hat <- lm(sample_eigenvalues[bulk_indices] ~ simulated_quantiles[bulk_indices] - 1)$coef[[1]]
  return(list(theta_hat = theta_hat, sigma2_hat = sigma2_hat))

}

#' @importFrom irlba irlba
#' @export
step2 <- function(sample_eigenvalues, n, p, beta, step1_output, M = 500,
                  use_irlba=TRUE)
{
  theta <- step1_output$theta_hat
  sigma2 <- step1_output$sigma2_hat

  l1_distribution <- sapply(1:M, function(i){
    Sigma_half <- diag(sqrt(sigma2 * rgamma(n = p, shape = theta, rate = theta)))
    X_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
    X <- X_raw %*% Sigma_half

    if (use_irlba)
    {
      if (n < p)
        1/n * irlba(X, nu=1, nv=0)$d^2
      else
        1/n * irlba(X, nu=0, nv=1)$d^2
    }
    else
    1/n * svd(X, nu=0, nv=0)$d[1]^2

  })

  quant <- quantile(l1_distribution, probs = 1 - beta)
  return(sum(sample_eigenvalues > quant))


}
