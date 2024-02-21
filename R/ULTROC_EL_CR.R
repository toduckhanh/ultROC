#' @export
ll_emp_2s <- function(n, p0, p_est, id_min_max){
  temp1 <- (p_est * log(p_est/p0) + ((1 - p_est) * log((1 - p_est)/(1 - p0))))
  temp1 <- 2 * c(n[1], n[id_min_max + 1]) * temp1
  if (any(is.nan(temp1))) res <- Inf
  else res <- sum(temp1)
  return(res)
}

#' @export
EL_cr_LTROC <- function(cpt, Xlist, xlim, ylim) {
  Sp_emp <- mean(Xlist[[1]] <= cpt)
  all_Se_emp <- sapply(Xlist[-1], function(y) {
    mean(y > cpt)
  })
  LSe_emp <- min(all_Se_emp)
  ### construct empirical likelihood statistics
  id_min_Se <- which.min(all_Se_emp)
  x_Sp <- seq(xlim[1], xlim[2], length.out = 300)
  y_Se <- seq(ylim[1], ylim[2], length.out = 300)
  n <- sapply(Xlist, length)
  ll_Sp_LSe <- sapply(x_Sp, function(x) {
    sapply(y_Se, function(y) {
      ll_emp_2s(n = n, p0 = c(x, y),
                p_est = c(1 - Sp_emp, LSe_emp),
                id_min_max = id_min_Se)
    })
  })
  return(list(Sp = Sp_emp, LSe = LSe_emp, x_Sp = x_Sp, y_Se = y_Se,
              ll_Sp_LSe = t(ll_Sp_LSe)))
}

#' @export
EL_cr_UTROC <- function(cpt, Xlist, xlim, ylim) {
  Sp_emp <- mean(Xlist[[1]] <= cpt)
  all_Se_emp <- sapply(Xlist[-1], function(y) {
    mean(y > cpt)
  })
  USe_emp <- max(all_Se_emp)
  ### construct empirical likelihood statistics
  id_max_Se <- which.max(all_Se_emp)
  x_Sp <- seq(xlim[1], xlim[2], length.out = 300)
  y_Se <- seq(ylim[1], ylim[2], length.out = 300)
  n <- sapply(Xlist, length)
  ll_Sp_USe <- sapply(x_Sp, function(x) {
    sapply(y_Se, function(y) {
      ll_emp_2s(n = n, p0 = c(x, y),
                p_est = c(1 - Sp_emp, USe_emp),
                id_min_max = id_max_Se)
    })
  })
  return(list(Sp = Sp_emp, USe = USe_emp, x_Sp = x_Sp, y_Se = y_Se,
              ll_Sp_USe = t(ll_Sp_USe)))
}
