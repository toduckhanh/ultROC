
#' @export
LTROC_emp <- function(p1, Xlist){
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmin, all_Se)
  return(res)
}

#' @export
LTROC_emp_2 <- function(Xlist, tt){
  Sp_emp <- mean(X[[1]] <= tt)
  Se_emp_all <- round(mapply(function(a) mean(a > tt), a = X[-1]), 7)
  id_min_Se <- which.min(Se_emp_all)
  LSe_emp <- min(Se_emp_all)
  return(c(Sp_emp, LSe_emp, id_min_Se))
}

#' @export
UTROC_emp <- function(p1, Xlist){
  t1 <- quantile(x = Xlist[[1]], probs = 1 - p1, names = FALSE)
  Fn <- lapply(Xlist[-1], function(x) ecdf(x))
  all_Se <- lapply(Fn, function(ff) 1 - ff(t1))
  res <- Reduce(pmax, all_Se)
  return(res)
}

#' @export
UTROC_emp_2 <- function(Xlist, tt){
  Sp_emp <- mean(X[[1]] <= tt)
  Se_emp_all <- round(mapply(function(a) mean(a > tt), a = X[-1]), 7)
  id_max_Se <- which.max(Se_emp_all)
  USe_emp <- max(Se_emp_all)
  return(c(Sp_emp, USe_emp, id_max_Se))
}

#' @export
ll_emp_2s <- function(n, p0, p_est, id_min_max){
  temp1 <- (p_est * log(p_est/p0) + ((1 - p_est) * log((1 - p_est)/(1 - p0))))
  temp1 <- 2 * c(n[1], n[id_min_max + 1]) * temp1
  if (any(is.nan(temp1))) res <- Inf
  else res <- sum(temp1)
  return(res)
}

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

#' @export
ll_emp_1s <- function(n, p0, p_est){
  temp1 <- (p_est * log(p_est/p0) + ((1 - p_est) * log((1 - p_est)/(1 - p0))))
  temp1 <- 2 * n * temp1
  if (is.nan(temp1)) res <- Inf
  else res <- temp1
  return(res)
}

#' @export
boot_LTAUC_emp <- function(X_list, p_est, B){
  ## resampling process
  empi_bts <- sapply(1:B, function(i){
    flag <- 0
    while(flag == 0){
      X_b <- lapply(X_list, function(x){
        sample(x, size = length(x), replace = TRUE)
      })
      mu_X_b <- sapply(X_b, mean)
      flag <- 1 - as.numeric(any(mu_X_b[1] > mu_X_b[-1]))
    }
    LTAUC_emp_bts <- LTAUC_emp(X_b)
    res <- ll_emp_1s(n = length(X_list[[1]]), p0 = p_est, p_est = LTAUC_emp_bts)
    return(res)
  })
  empi_bts[is.na(empi_bts)] <- Inf
  r_est <- ((7 / 9)^3) / median(empi_bts)
  return(r_est)
}

#' @export
boot_UTAUC_emp <- function(X_list, p_est, B){
  ## resampling process
  empi_bts <- sapply(1:B, function(i){
    flag <- 0
    while(flag == 0){
      X_b <- lapply(X_list, function(x){
        sample(x, size = length(x), replace = TRUE)
      })
      mu_X_b <- sapply(X_b, mean)
      flag <- 1 - as.numeric(any(mu_X_b[1] > mu_X_b[-1]))
    }
    UTAUC_emp_bts <- UTAUC_emp(X_b)
    res <- ll_emp_1s(n = length(X_list[[1]]), p0 = p_est, p_est = UTAUC_emp_bts)
    return(res)
  })
  empi_bts[is.na(empi_bts)] <- Inf
  r_est <- ((7 / 9)^3) / median(empi_bts)
  return(r_est)
}

#' @export
bts_func <- function(Y, n, LTAUC_est, UTAUC_est, seed, B = 200){
  if (missing(seed)) seed <- 34
  set.seed(seed)
  res_bts <- sapply(1:B, function(i){
    flag <- 0
    while(flag == 0){
      Y_b <- lapply(Y, function(x){
        sample(x, size = length(x), replace = TRUE)
      })
      mu_Y <- sapply(Y_b, mean)
      flag <- all(mu_Y[1] < mu_Y[-1])
    }
    LTAUC_bts <- LTAUC_emp(Y_b)
    UTAUC_bts <- UTAUC_emp(Y_b)
    if (LTAUC_bts == 1) LTAUC_bts <- n[1] / (n[1] + 0.5)
    if (UTAUC_bts == 1) UTAUC_bts <- n[1] / (n[1] + 0.5)
    # SIK_bts <- SIK_hist(Y_b)
    res_1 <- ll_emp_1s(n = n[1], p0 = LTAUC_est, p_est = LTAUC_bts)
    res_2 <- ll_emp_1s(n = n[1], p0 = UTAUC_est, p_est = UTAUC_bts)
    # res_5 <- SIK_bts
    res_3 <- LTAUC_bts
    res_4 <- UTAUC_bts
    return(c(res_1, res_2, res_3, res_4))
  })
  return(res_bts)
}

myfun <- function(theta, theta_est, r_adj, qc, n) {
  ll_est <- ll_emp_1s(n = n, p0 = theta, p_est = theta_est)
  if (is.na(ll_est)) ll_est <- Inf
  ll_est_adj <- ll_est
  if (!is.infinite(ll_est)){
    ll_est_adj <- r_adj * ll_est_adj
  }
  return(ll_est_adj - qc)
}

#' @export
EL_ci_LUTAUC <- function(LUTAUC_est, n, w_est, ci_level){
  res <- lapply(ci_level, function(x){
    LI_eng <- uniroot(f = myfun, interval = c(0, LUTAUC_est),
                      theta_est = LUTAUC_est, qc = qchisq(x, 1), r_adj = w_est,
                      n = n)$root
    UI_eng <- uniroot(f = myfun, interval = c(LUTAUC_est, 1),
                      theta_est = LUTAUC_est, qc = qchisq(x, 1), r_adj = w_est,
                      n = n)$root
    ci_vus <- c(LI_eng, UI_eng)
  })
  return(res)
}
