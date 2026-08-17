
#' Fit the REMARCS model with a mixed-effects regression component
#'
#' Implements the REMARCS algorithm using an initial mixed-effects model
#' and a latent ARMA(1, 1) time-series component.
#'
#' @param mod The initial mixed-effects model fitted with `lme4::lmer()`.
#' @param dat The data used for fitting. The variable for time must be named
#'   `t`.
#' @param resp_var Character string giving the name of the response variable.
#' @param tol Convergence tolerance for the iterative REMARCS algorithm.
#'
#' @return A list containing the estimated regression parameters, time-series
#'   parameters, log-likelihood, AIC, and fitted values.
#'
#' @export
RCSmixed.fit <- function(mod, dat, resp_var, tol = 0.1, max_updates = 10, arma_order = c(1,1)) {




  # --------------------------------------------------------------------------
  # Initial mixed-effects model and daily residual summaries
  # --------------------------------------------------------------------------
  y = dat[[resp_var]]

  resid <- stats::residuals(mod)
  dat$resid <- resid
  var_re = sigma(mod)^2 * tcrossprod(getME(mod, "Lambda"))

  # --------------------------------------------------------------------------
  # Aggregate residuals
  # --------------------------------------------------------------------------
  df_w <- dat |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      w = mean(resid),
      S2 = stats::var(resid),
      n = dplyr::n(),
      .groups = "drop"
    )
  n_t = df_w$n


  # Fill in missing days so that the time-series component has a complete
  # sequence of time points.

  df_w <- dplyr::left_join(
    tibble::tibble(t = 1:max(df_w$t)),
    df_w,
    by = "t"
  )

  stats::acf(
    df_w$w,
    na.action = stats::na.pass,
    main = "ACF of average residuals"
  )

  # --------------------------------------------------------------------------
  # Initial ARMA(1, 1) estimates
  # --------------------------------------------------------------------------

  arma_w <- stats::arima(
    df_w$w,
    order = c(1, 0, 1),
    include.mean = FALSE
  )

  print(acf(df_w$w, na.action = na.pass))
  print(arma_w)
  # S_o = Sigma_o_varyingn(ar = phi_hat, ma_b = theta_b_hat, veps = var_eps_hat, v_b = var_b_hat, n_t = n_t)
  # conf_region_2d(theta_hat = c(phi_hat, theta_b_hat),
  #                theta_true = c(phi_hat, theta_b_hat),
  #                Sigma = S_o[1:2,1:2] / T) |> print()


  phi_hat <- arma_w$coef["ar1"]
  theta_b_hat <- arma_w$coef["ma1"]
  var_b_hat <- arma_w$sigma2
  ans1 = readline(prompt = "Would you like to restrict the ARMA model?
                 [1] AR(1)
                 [2] MA(1)
                 [3] No
                 ")

  if (ans1 == 1){
    arma_w <- stats::arima(
      df_w$w,
      order = c(1, 0, 0),
      include.mean = FALSE
    )
    phi_hat = arma_w$coef["ar1"]
    theta_b_hat = 0
    var_eta_hat = arma_w$sigma2
  }else if(ans1 == 2){
    arma_w <- stats::arima(
      df_w$w,
      order = c(0, 0, 1),
      include.mean = FALSE
    )
    phi_hat = 0
    theta_b_hat = arma_w$coef["ma1"]
    var_eta_hat = arma_w$sigma2
  }



  df_w$n <- tidyr::replace_na(df_w$n, 0)
  df_w_initial <- df_w

  # --------------------------------------------------------------------------
  # Design matrices and initial variance estiamtes
  # --------------------------------------------------------------------------

  X_f <- stats::model.matrix(mod)
  X_f_m <- X_f

  X_r <- lme4::getME(mod, "Z")
  vc <- lme4::VarCorr(mod)

  var_re <- stats::sigma(mod)^2 *
    Matrix::tcrossprod(lme4::getME(mod, "Lambda"))

  l1 <- as.numeric(stats::logLik(mod))
  dat$resid_initial <- dat$resid

  n_t <- df_w$n
  X_t <- f_Xt(n_t)
  V <- f_V(n_t)
  T <- length(n_t)
  N <- diag(n_t)

  df_w_veps_est <- df_w |>
    dplyr::filter(n > 1)

  var_eps_hat <- sum(
    (df_w_veps_est$n - 1) * df_w_veps_est$S2
  ) / (
    sum(df_w_veps_est$n) - length(df_w_veps_est$n)
  )

  var_a_hat <- var_eps_hat * mean((n_t[n_t != 0])^(-1))

  # --------------------------------------------------------------------------
  # Apply ARMAN method for latent parameter point estimates
  # --------------------------------------------------------------------------

  wong_est <- wong(
    ar_hat = phi_hat,
    ma_b_hat = theta_b_hat,
    v_b_hat = var_b_hat,
    v_a_hat = var_a_hat,
    theta_free = arma_order[2] == 1
  )

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]


  # --------------------------------------------------------------------------
  # Estimate Time Series Covariance Matrix
  # --------------------------------------------------------------------------

  G_hat <- fG(
    size = T,
    eta_variance = var_eta_hat,
    ar = phi_hat,
    ma = theta_hat
  )

  # --------------------------------------------------------------------------
  # Estimate Time Series and calculate likelihood
  # --------------------------------------------------------------------------

  u_hat <- f_u_hat(
    n_t,
    df_w$w,
    G_hat,
    var_eps_hat
  )

  observed_days <- n_t != 0
  G_obs <- G_hat[observed_days, observed_days]

  l2 <- -1 / 2 * (
    T * log(2 * pi) +
      as.numeric(
        determinant(
          G_obs,
          logarithm = TRUE
        )$modulus
      ) +
      as.vector(
        t(u_hat) %*%
          solve(G_obs) %*%
          u_hat
      )
  )

  # Current likelihood
  l_new <- l2 + l1
  l_old <- l_new + 2 * tol

  # --------------------------------------------------------------------------
  # Set up vectors to track parameter estimates and likelihoods
  # --------------------------------------------------------------------------

  l1_track <- l1
  l2_track <- l2
  l_track <- l_old

  phi_track <- phi_hat
  theta_track <- theta_hat
  var_eta_track <- var_eta_hat
  var_eps_track <- var_eps_hat



  # --------------------------------------------------------------------------
  # REMARCS optimization loop
  # --------------------------------------------------------------------------
  updates = 1
  while (abs(l_old - l_new) > tol & updates <= max_updates) {

    # --------------------------------------------------------------------------
    # track updates
    # --------------------------------------------------------------------------
    updates = updates + 1
    l_old <- l_new


    # --------------------------------------------------------------------------
    # Generate transformed eigendecomposition to facilitate later calculations
    # NOTE: This may prove later to be unnecessary
    # --------------------------------------------------------------------------
    NGN_eigen <- eigen(
      N^(1 / 2) %*% G_hat %*% N^(1 / 2)
    )

    NGN_evec <- NGN_eigen$vectors
    NGN_eval <- NGN_eigen$values
    L <- diag(NGN_eval)

    V <- f_V(n_t)

    # --------------------------------------------------------------------------
    # Extract time series estimate and refit mixed effect model to update conditional likelihood
    # Eventually unnecessary if conditioning on aggregate values instead of time series
    # --------------------------------------------------------------------------
    dat$less_u_hat <-
      dat[[resp_var]] -
      as.vector(X_t[, observed_days] %*% u_hat)

    new_formula = formula(mod)
    new_formula[[2]] = quote(less_u_hat)
    re_new = lme4::lmer(data = dat, formula = new_formula)

    l1 <- as.numeric(stats::logLik(re_new))



    #Sparse versions of the large design matrices.
    X_t_sparse <- Matrix::Matrix(X_t, sparse = TRUE)
    X_r_sparse <- Matrix::Matrix(X_r, sparse = TRUE)

    # ------------------------------------------------------------------------
    # Woodbury calculations are used here to lessen calculation time
    # ------------------------------------------------------------------------

    G_inv <- chol2inv(chol(G_hat))

    Gplus_inv <- chol2inv(
      chol(G_inv + N / var_eps_hat)
    )

    B_inv <- (solve(var_re) + Matrix::t(X_r) %*% X_r / var_eps_hat - var_eps_hat ^ (-2) * Matrix::t(X_r) %*% X_t_sparse %*% Gplus_inv %*% Matrix::t(X_t_sparse) %*% X_r_sparse) #|>
    B_inv <- chol2inv(chol(B_inv))

    beta_hat_2 = estimate_beta_gls(y = y,
                                   X_f = X_f_m,
                                   X_r = X_r_sparse,
                                   X_t = X_t_sparse,
                                   Sigma_b = var_re,
                                   Gplus_inv = Gplus_inv,
                                   sigma_eps2 = var_eps_hat)
    beta_hat = beta_hat_2$beta_hat
    # ------------------------------------------------------------------------
    # Update random effect intercepts and variance estimates
    # ------------------------------------------------------------------------

    b_hat_list = estimate_b(y = y,
                            X_f = X_f_m,
                            X_r = X_r_sparse,
                            X_t = X_t_sparse,
                            beta_hat = beta_hat,
                            Sigma_b = var_re,
                            Gplus_inv = Gplus_inv,
                            sigma_eps2 = var_eps_hat)

    b_hat <- b_hat_list$b_hat
    B <- b_hat_list$B

    B_inv <- solve(B)
    var_re_update <- update_Sigma_b(b_hat = b_hat, B_inv = B_inv, mod = mod)


    var_re <- var_re_update$matrix



    # ------------------------------------------------------------------------
    # Update ARMA-Noise Process Estimate
    # ------------------------------------------------------------------------


    dat$resid <-
      as.vector(
        dat[[resp_var]] -
          X_f_m %*% beta_hat - X_r %*% b_hat
      )

    # ------------------------------------------------------------------------
    # Repeat fitting of ARMA-Noise process
    # ------------------------------------------------------------------------


    df_w <- dat |>
      dplyr::group_by(t) |>
      dplyr::summarise(
        w = mean(resid),
        S2 = stats::var(resid),
        n = dplyr::n(),
        .groups = "drop"
      )

    df_w <- dplyr::left_join(
      tibble::tibble(t = seq_len(max(df_w$t))),
      df_w,
      by = "t"
    )

    df_w$n <- tidyr::replace_na(df_w$n, 0)

    arma_w <- forecast::Arima(
      df_w$w,
      order = c(1, 0, 1),
      include.mean = FALSE,
      #init = c(phi_hat, theta_b_hat),
      method = "ML"
    )

    phi_hat <- arma_w$coef["ar1"]
    theta_b_hat <- arma_w$coef["ma1"]
    var_b_hat <- arma_w$sigma2



    df_w_veps_est <- df_w |>
      dplyr::filter(n > 1)

    var_eps_hat <- sum(
      (df_w_veps_est$n - 1) * df_w_veps_est$S2
    ) / (
      sum(df_w_veps_est$n) - length(df_w_veps_est$n)
    )

    var_a_hat <- var_eps_hat * mean((n_t[n_t != 0])^(-1))

    wong_est <- wong(
      phi_hat,
      theta_b_hat,
      var_b_hat,
      var_a_hat,
      theta_free = arma_order[2] == 1
    )

    theta_hat <- wong_est[1]
    var_eta_hat <- wong_est[2]

    # ------------------------------------------------------------------------
    # Update estimates of time series and corresponding covariance matrix
    # ------------------------------------------------------------------------

    G_hat <- fG(
      size = T,
      eta_variance = var_eta_hat,
      ar = phi_hat,
      ma = theta_hat
    )

    u_hat <- f_u_hat(
      n_t,
      df_w$w,
      G_hat,
      var_eps_hat
    )

    G_obs <- G_hat[observed_days, observed_days]

    # ------------------------------------------------------------------------
    # Update level 2 likelihood
    # ------------------------------------------------------------------------


    l2 <- -1 / 2 * (
      T * log(2 * pi) +
        as.numeric(
          determinant(
            G_obs,
            logarithm = TRUE
          )$modulus
        ) +
        as.vector(
          t(u_hat) %*%
            solve(G_obs) %*%
            u_hat
        )
    )

    # cat("l2 = ", l2, "\n")

    l_new <- l2 + l1

    # ------------------------------------------------------------------------
    # Store parameter estimates and likelihoods
    # ------------------------------------------------------------------------


    l1_track <- append(l1_track, l1)
    l2_track <- append(l2_track, l2)
    l_track <- append(l_track, l_new)


    phi_track <- append(phi_track, phi_hat)
    theta_track <- append(theta_track, theta_hat)
    var_eta_track <- append(var_eta_track, var_eta_hat)
    var_eps_track <- append(var_eps_track, var_eps_hat)

  }

  # --------------------------------------------------------------------------
  # Asymptotic covariance of time-series parameter estimates
  # --------------------------------------------------------------------------

  S_o <- Sigma_o_varyingn(
    phi_hat,
    theta_b_hat,
    var_eps_hat,
    var_b_hat,
    n_t
  )

  wong_est <- wong(
    phi_hat,
    theta_b_hat,
    var_b_hat,
    var_eps_hat / mean(n_t[n_t != 0]),
    theta_free = arma_order[2] == 1
  )

  mean_inv_n = mean(n_t[n_t != 0])

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]
  var_eta_est <- var_eta_hat


  if (arma_order[2] == 1){
  H1 <- dHf_dtau1(
    theta_hat,
    var_eta_hat
  )

  H0 <- dHf_dtau0(
    phi_hat,
    theta_b_hat,
    var_a_hat,
    var_b_hat
  )

  GGG <- -solve(H1) %*% H0
  S_f <- GGG %*% S_o %*% Matrix::t(GGG)
  S_l <- S_f[c(1, 5, 6), c(1, 5, 6)]
  S_l_plot <- S_l / T

  ts.pars <- c(
    phi_hat,
    theta_hat
  )

  names(ts.pars) <- c(
    "phi",
    "theta"
  )

  ts.pars <- cbind(
    ts.pars,
    sqrt(diag(S_l_plot)[1:2])
  )

  pvals_ts <- 2 * (
    1 - stats::pnorm(
      abs(ts.pars[, 1] / ts.pars[, 2])
    )
  )

  ts.pars <- cbind(
    ts.pars,
    pvals_ts
  )

  print(ts.pars)

  } else{
  #var_eta_est = max((1 + theta_b_hat ^ 2) * var_b_hat - (1 + phi_hat ^ 2) * var_a_hat, 0)
  GG = matrix(0, nrow = 3, ncol = 4)
  GG[1,1] = 1
  GG[2,4] = 1
  GG[3,1] = -2 * phi_hat * var_a_hat
  GG[3,2] = 2 * theta_b_hat * var_b_hat
  GG[3,3] = 1 + theta_b_hat ^ 2
  GG[3,4] = -(1 + phi_hat ^ 2)
  S_f = GG %*% S_o %*% t(GG)

  S_l_plot = S_f / T
  #var_eta_hat = var_eta_est
  #theta_hat = 0

  ts.pars <- c(
    phi_hat
  )

  names(ts.pars) <- c(
    "phi"
  )

  ts.pars <- cbind(
    ts.pars,
    sqrt(diag(S_l_plot)[c(1)])
  )

  pvals_ts <- 2 * (
    1 - stats::pnorm(
      abs(ts.pars[, 1] / ts.pars[, 2])
    )
  )

  ts.pars <- cbind(
    ts.pars,
    pvals_ts
  )

  print(ts.pars)
  }








  Sigma_beta = var_beta(X_f,
                        X_r,
                        X_t,
                        Sigma_b = var_re,
                        Gplus_inv = Gplus_inv,
                        sigma_eps2 = var_eps_hat)




  # --------------------------------------------------------------------------
  # Final estimates and fitted values
  # --------------------------------------------------------------------------

  beta_hat <- cbind(
    beta_hat,
    sqrt(diag(Sigma_beta))
  )

  pvals <- 2 * (
    1 - stats::pnorm(
      abs(beta_hat[, 1] / beta_hat[, 2])
    )
  )

  beta_hat <- cbind(
    beta_hat,
    pvals
  )

  AIC <- 2 * (nrow(beta_hat) + 3) - 2 * l_new

  fitted.mat <- data.frame(
    Day = dat$t,
    Reg = X_f %*% beta_hat[, 1],
    TS = X_t[, observed_days] %*% u_hat
  )


  m1 <- nlevels(getME(mod, "flist")$adm1)
  m2 <- nlevels(getME(mod, "flist")$`adm1:adm2`)

  ranef_final = ranef(mod)
  ranef_final$adm1 = b_hat[(m2 + 1):(m1 + m2)]
  ranef_final$adm1 = as.matrix(ranef_final$adm1)
  ranef_final$`adm1:adm2` = b_hat[(1):(m2)]
  ranef_final$`adm1:adm2` = as.matrix(ranef_final$`adm1:adm2`)

  rownames(ranef_final$adm1) = rownames(ranef(mod)$adm1)
  rownames(ranef_final$`adm1:adm2`) = rownames(ranef(mod)$`adm1:adm2`)

  res.list <- list(
    beta_hat,
    b_hat = ranef_final,
    var_b = var_re_update$variances,
    ts.pars,
    l1_track,
    l2_track,
    AIC,
    fitted.mat,
    var_eta_hat
  )

dat$final_resid = dat[[resp_var]] - fitted.mat[,2] - fitted.mat[,3]
final_resid_daily = dat |>
  group_by(t) |>
  summarize(mean_resid = mean(final_resid))

final_resid_daily = left_join(tibble(t = 1:max(final_resid_daily$t)), final_resid_daily)

print(acf(final_resid_daily$mean_resid,na.action = na.pass, main = "ACF of Final Aggregated Residuals Aggregated"))


  names(res.list) <- c(
    "Beta_Pars",
    "Random Effects",
    "Random Effect Variances",
    "TS_Pars",
    "level 1 log(likelihood)",
    "level 2 log(likelihood)",
    "AIC",
    "Fitted_Values",
    "Time Series Variance"
  )
  ans2 = readline("Would you like to refit with a latent AR(1) process?
                  [1] Yes
                  [2] No")

  if (ans2 == 1){
    return(RCSmixed.fit(mod = mod,
                        dat = dat,
                        resp_var = resp_var,
                        tol = tol,
                        max_updates = max_updates,
                        arma_order = c(1,0)))
  }else{
  return(res.list)
  }
}

