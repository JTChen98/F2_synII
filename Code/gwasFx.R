gwas <- function (y, X)
{
  n <- length(y)
  U1 = sum(y)
  U2 = U1/n
  ytr = y - rep(1, n) * U2
  U3 = colSums(X)
  U4 = U3/n
  Str = X - tcrossprod(rep(1, n), U4)
  Str2 = colSums(Str^2)
  b = as.vector(crossprod(ytr, Str)/Str2)
  sig = (sum(ytr^2) - b^2 * Str2)/(n - 2)
  err = sqrt(sig * (1/Str2))
  pval = pchisq(b^2/err^2, 1, low = FALSE)
  return(data.frame(beta = b, se = err, p = pval))
}
