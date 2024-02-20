#' @export
LTROC_emp <- function(p1, Xlist){
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmin, all_Se)
  return(res)
}

#' @export
UTROC_emp <- function(p1, Xlist){
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmax, all_Se)
  return(res)
}
