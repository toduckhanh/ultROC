#' @import stats

#' @export
ll_emp_1s <- function(n, p0, p_est){
  temp1 <- (p_est * log(p_est/p0) + ((1 - p_est) * log((1 - p_est)/(1 - p0))))
  temp1 <- 2 * n * temp1
  if (is.nan(temp1)) res <- Inf
  else res <- temp1
  return(res)
}

#' @export
bts_func <- function(Y, n, LTAUC_est, UTAUC_est, seed, B = 200){
  if (missing(seed)) seed <- 34
  set.seed(seed)
  res_bts <- sapply(1:B, function(i){
    # flag <- 0
    # while(flag == 0){
      Y_b <- lapply(Y, function(x){
        sample(x, size = length(x), replace = TRUE)
      })
      # mu_Y <- sapply(Y_b, mean)
      # flag <- all(mu_Y[1] < mu_Y[-1])
    # }
    LTAUC_bts <- LTAUC_emp(Y_b)
    UTAUC_bts <- UTAUC_emp(Y_b)
    if (LTAUC_bts == 1) LTAUC_bts <- n[1] / (n[1] + 0.5)
    if (UTAUC_bts == 1) UTAUC_bts <- n[1] / (n[1] + 0.5)
    res_1 <- ll_emp_1s(n = n[1], p0 = LTAUC_est, p_est = LTAUC_bts)
    res_2 <- ll_emp_1s(n = n[1], p0 = UTAUC_est, p_est = UTAUC_bts)
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
