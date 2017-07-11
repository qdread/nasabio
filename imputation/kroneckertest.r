# Test Kronecker product code

Sigma <- matrix(rnorm(9),nrow=3)
Sigma <- Sigma %*% t(Sigma)
R <- diag(2)

m <- nrow(Sigma)
n <- nrow(R)

lambda_correct <- Sigma %x% R

lambda <- matrix(0, nrow = 6, ncol = 6)

for (i in 1:m)
  for (j in 1:m)
    for (k in 1:n)
      for (l in 1:n)
        lambda[n*(i-1)+k, n*(j-1)+l] = Sigma[i, j] * R[k, l];