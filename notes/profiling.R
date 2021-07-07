library(microbenchmark)

n <- 500
p <- 100
K <- 5
rho <- 5
mu_k <- rho * sqrt(p/n)
vectors <- svd(matrix(rnorm(p*p), p, p))$u[,1:K]
theta <- 5
sigma2 <- 1
Sigma <- tcrossprod(mu_k * vectors) + diag(rgamma(p, theta, rate = theta/sigma2))

X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
S <- 1/n * tcrossprod(X)
step1_out <- step1(eigen(S)$values, n, p, 0.1, B = 20)
step2(eigen(S)$values, n, p, beta = 0.2, M = 500, step1_output = step1_out)
