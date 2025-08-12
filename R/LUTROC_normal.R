
#' @export
LTROC_normal <- function(p1, mu, sigma){
  t1 <- qnorm(1 - p1, mean = mu[1], sd = sigma[1])
  all_Se <- mapply(function(aa, bb) {
    pnorm(t1, mean = aa, sd = bb, lower.tail = FALSE)
  }, aa = mu[-1], bb = sigma[-1], SIMPLIFY = FALSE)
  res <- Reduce(pmin, all_Se)
  return(res)
}

#' @export
UTROC_normal <- function(p1, mu, sigma){
  t1 <- qnorm(1 - p1, mean = mu[1], sd = sigma[1])
  all_Se <- mapply(function(aa, bb) {
    pnorm(t1, mean = aa, sd = bb, lower.tail = FALSE)
  }, aa = mu[-1], bb = sigma[-1], SIMPLIFY = FALSE)
  res <- Reduce(pmax, all_Se)
  return(res)
}

#' @export
LTAUC_normal <- function(mu, sigma){
  integrate(LTROC_normal, lower = 0, upper = 1, mu = mu, sigma = sigma,
            rel.tol = 1e-9)$value
}

#' @export
UTAUC_normal <- function(mu, sigma){
  integrate(UTROC_normal, lower = 0, upper = 1, mu = mu, sigma = sigma,
            rel.tol = 1e-9)$value
}
