#' @export
LTROC_emp <- function(Xlist, p1){
  if (missing(p1)) p1 <- seq(0, 1, by = 0.001)
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmin, all_Se)
  res[p1 == 0] <- 0
  res[p1 == 1] <- 1
  return(cbind(fpr = p1, lse = res))
}

#' @export
UTROC_emp <- function(Xlist, p1){
  if (missing(p1)) p1 <- seq(0, 1, by = 0.001)
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmax, all_Se)
  res[p1 == 0] <- 0
  res[p1 == 1] <- 1
  return(cbind(fpr = p1, use = res))
}
