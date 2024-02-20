#' @import utils

#' @export
LTAUC_emp <- function(X_list){
  ## X_list: list of X1, X2, ..., Xk
  xx <- X_list[[1]]
  Fn <- lapply(X_list[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(xx))
  LSe_est <- Reduce(pmin, all_Se)
  return(mean(LSe_est))
}

#' @export
UTAUC_emp <- function(X_list){
  ## X_list: list of X1, X2, ..., Xk
  xx <- X_list[[1]]
  Fn <- lapply(X_list[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(xx))
  USe_est <- Reduce(pmax, all_Se)
  return(mean(USe_est))
}
