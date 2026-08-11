
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
RCSmixed.fit <- function(mod, dat, resp_var, tol = 0.1, max_updates = 10) {

  #### Cheating with fixed.fit
  # rand_eff <- predict(mod, re.form = NULL) -
  #   predict(mod, re.form = NA)
  #
  # y_star <- model.response(model.frame(mod)) - rand_eff
  #
  # fixed_form = lme4::nobars(formula(mod))
  #
  # dat$less_re = y_star
  #
  # fixed_form[[2]] = quote(less_re)
  # no_re = lm(data = dat, formula = fixed_form)
  # fixed.fit = remarcs::RCSfixed.fit(mod = no_re, dat = dat, resp_var = "less_re", tol)
  # return(fixed.fit)

  ######

  # --------------------------------------------------------------------------
  # Initial mixed-effects model and daily residual summaries
  # --------------------------------------------------------------------------
  y = dat[[resp_var]]

  resid <- stats::residuals(mod)
  dat$resid <- resid
  var_re = sigma(mod)^2 * tcrossprod(getME(mod, "Lambda"))

  df_w <- dat |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      w = mean(resid),
      S2 = stats::var(resid),
      n = dplyr::n(),
      .groups = "drop"
    )
  n_t = df_w$n

  # Plot the original response. This is retained from the original function.
  # plot.orig <- ggplot2::ggplot(data = dat) +
  #   ggplot2::geom_point(
  #     mapping = ggplot2::aes(x = t, y = dat[[resp_var]])
  #   )

  # y <- dat[[resp_var]]

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

  phi_hat <- arma_w$coef["ar1"]
  theta_b_hat <- arma_w$coef["ma1"]
  var_b_hat <- arma_w$sigma2



  df_w$n <- tidyr::replace_na(df_w$n, 0)
  df_w_initial <- df_w

  # --------------------------------------------------------------------------
  # Design matrices and initial mixed-effects quantities
  # --------------------------------------------------------------------------

  X_f <- stats::model.matrix(mod)
  X_f_m <- X_f

  X_r <- lme4::getME(mod, "Z")
  vc <- lme4::VarCorr(mod)

  var_re <- stats::sigma(mod)^2 *
    Matrix::tcrossprod(lme4::getME(mod, "Lambda"))

  l1 <- as.numeric(stats::logLik(mod))
  dat$resid_initial <- dat$resid

  # --------------------------------------------------------------------------
  # Initial variance estimates
  # --------------------------------------------------------------------------

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

  wong_est <- wong(
    ar_hat = phi_hat,
    ma_b_hat = theta_b_hat,
    v_b_hat = var_b_hat,
    v_a_hat = var_a_hat
  )

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]

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
    updates = updates + 1
    l_old <- l_new

    NGN_eigen <- eigen(
      N^(1 / 2) %*% G_hat %*% N^(1 / 2)
    )

    NGN_evec <- NGN_eigen$vectors
    NGN_eval <- NGN_eigen$values
    L <- diag(NGN_eval)

    V <- f_V(n_t)

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
    # Woodbury calculations
    # ------------------------------------------------------------------------

    G_inv <- chol2inv(chol(G_hat))

    Gplus_inv <- chol2inv(
      chol(G_inv + N / var_eps_hat)
    )

    B_inv <- (solve(var_re) + Matrix::t(X_r) %*% X_r / var_eps_hat - var_eps_hat ^ (-2) * Matrix::t(X_r) %*% X_t_sparse %*% Gplus_inv %*% Matrix::t(X_t_sparse) %*% X_r_sparse) #|>
    B_inv <- chol2inv(chol(B_inv))


    R1 <-
      Matrix::crossprod(X_f_m) / var_eps_hat -
      var_eps_hat^(-2) *
      Matrix::crossprod(X_f_m, X_t) %*%
      Gplus_inv %*%
      Matrix::crossprod(X_t, X_f_m)

    Bchol <- t(chol(B_inv))

    R2_half <-
      Matrix::t(X_f_m) / var_eps_hat -
      var_eps_hat^(-2) *
      Matrix::t(X_f_m) %*%
      X_t_sparse %*%
      Gplus_inv %*%
      Matrix::t(X_t_sparse)

    R2_half <- R2_half %*% X_r_sparse %*% Bchol
    R2 <- R2_half %*% Matrix::t(R2_half)

    R <- R1 + R2
    R_inv <- chol2inv(chol(R))

    # ------------------------------------------------------------------------
    # GLS estimate of beta
    # ------------------------------------------------------------------------

    S1 <-
      y / var_eps_hat -
      as.vector(
        var_eps_hat^(-2) *
          X_t_sparse %*%
          (
            Gplus_inv %*%
              (Matrix::t(X_t_sparse) %*% y)
          )
      )

    S2 <- Matrix::t(X_r_sparse) %*% S1
    S2 <- B_inv %*% S2
    S2 <- X_r_sparse %*% S2

    S2 <-
      S2 / var_eps_hat -
      var_eps_hat^(-2) *
      (
        X_t_sparse %*%
          (
            Gplus_inv %*%
              (Matrix::t(X_t_sparse) %*% S2)
          )
      )

    S <- S1 - S2

    beta_hat <- R_inv %*% Matrix::t(X_f_m) %*% S
    beta_hat_2 = estimate_beta_gls(y = y,
                                   X_f = X_f_m,
                                   X_r = X_r_sparse,
                                   X_t = X_t_sparse,
                                   Sigma_b = var_re,
                                   Gplus_inv = Gplus_inv,
                                   sigma_eps2 = var_eps_hat)
    beta_hat = beta_hat_2$beta_hat
    # ------------------------------------------------------------------------
    # Update residuals and time-series component
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


    var_re <- var_re_update$var_re





    dat$resid <-
      as.vector(
        dat[[resp_var]] -
          X_f_m %*% beta_hat - X_r %*% b_hat
      )

    ###DONT FORGET TO ADD UPDATE to Sigma_b######


    # dat$resid <- stats::residuals(re_new)
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
      var_a_hat
    )

    theta_hat <- wong_est[1]
    var_eta_hat <- wong_est[2]


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
    var_eps_hat / mean(n_t[n_t != 0])
  )

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]
  var_eta_est <- var_eta_hat

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

  vec_est <- c(
    phi_hat,
    theta_hat,
    var_eta_hat
  )

  vec_true <- c(
    phi_hat,
    theta_hat,
    var_eta_hat
  )

  vec_name <- c(
    "phi",
    "theta",
    "var_eta"
  )

  # --------------------------------------------------------------------------
  # Standard errors for beta
  # --------------------------------------------------------------------------

  XftXf <- Matrix::crossprod(X_f_m)
  XftXf_inv <- solve(XftXf)

  Sigma_beta3_term1 <-
    Matrix::t(X_f_m) %*%
    X_t[, observed_days] %*%
    diag(sqrt(1 / n_t[observed_days])) %*%
    NGN_evec[observed_days, observed_days] %*%
    diag(
      sqrt(
        1 / (
          NGN_eval[observed_days] +
            var_eps_hat
        )
      )
    )

  Sigma_beta3_term1 <-
    Sigma_beta3_term1 %*%
    Matrix::t(Sigma_beta3_term1)

  Sigma_beta3_term2 <-
    Matrix::t(X_f_m) %*%
    V

  Sigma_beta3_term2 <-
    Sigma_beta3_term2 %*%
    Matrix::t(Sigma_beta3_term2) /
    var_eps_hat

  Sigma_beta3 <-
    solve(
      Sigma_beta3_term1 +
        Sigma_beta3_term2
    )

  # --------------------------------------------------------------------------
  # Final estimates and fitted values
  # --------------------------------------------------------------------------

  beta_hat <- cbind(
    beta_hat,
    sqrt(diag(Sigma_beta3))
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

  ts.pars <- c(
    phi_hat,
    theta_hat,
    var_eta_hat
  )

  names(ts.pars) <- c(
    "phi",
    "theta",
    "var_eta"
  )

  res.list <- list(
    beta_hat,
    ts.pars,
    l_new,
    AIC,
    fitted.mat
  )

  names(res.list) <- c(
    "Beta_Pars",
    "TS_Pars",
    "log(likelihood)",
    "AIC",
    "Fitted_Values"
  )

  return(res.list)
}

