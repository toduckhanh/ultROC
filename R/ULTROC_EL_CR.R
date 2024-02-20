#' @export
ll_emp_2s <- function(n, p0, p_est, id_min_max){
  temp1 <- (p_est * log(p_est/p0) + ((1 - p_est) * log((1 - p_est)/(1 - p0))))
  temp1 <- 2 * c(n[1], n[id_min_max + 1]) * temp1
  if (any(is.nan(temp1))) res <- Inf
  else res <- sum(temp1)
  return(res)
}


