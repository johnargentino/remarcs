#' Fit the REMARCS model with a fixed- or mixed-effects regression component
#'
#' Implements the REMARCS algorithm using either an initial fixed-effects
#' model fitted with `lm()` or a mixed-effects model fitted with
#' `lme4::lmer()`, together with a latent ARMA time-series component.
#'
#' @param mod The initial regression model. Must be an `lm` or `merMod`
#'   model object.
#' @param dat The data used for fitting. The variable for time must be
#'   named `t`.
#' @param resp_var Character string giving the name of the response variable.
#' @param tol Convergence tolerance for the iterative REMARCS algorithm.
#' @param max_updates Maximum number of REMARCS updates.
#' @param arma_order Integer vector giving the AR and MA orders for the
#'   latent process. Currently `c(1, 1)` and `c(1, 0)` are supported.
#'
#' @return A list containing estimated regression parameters, time-series
#'   parameters, likelihood information, AIC, fitted values, and, for
#'   mixed models, random effects and their variance estimates.
#'
#' @export
RCS.fit <- function(mod,
                    dat,
                    resp_var,
                    tol = 0.1,
                    max_updates = 10,
                    arma_order = c(1, 1)) {

  # ==========================================================================

  # 1. Validate inputs and determine model type

  # ==========================================================================

  is_mixed <- inherits(mod, "merMod")
  is_fixed <- inherits(mod, "lm") && !is_mixed

  if (!is_fixed && !is_mixed) {
    stop(
      "`mod` must be either an `lm` model or an lme4 mixed-effects ",
      "model inheriting from `merMod`."
    )
  }

  if (!is.character(resp_var) || length(resp_var) != 1) {
    stop("`resp_var` must be a single character string.")
  }

  if (!resp_var %in% names(dat)) {
    stop("`resp_var` was not found in `dat`.")
  }

  if (!"t" %in% names(dat)) {
    stop("The time variable in `dat` must be named `t`.")
  }

  if (length(arma_order) != 2 ||
      !all(arma_order %in% c(0, 1)) ||
      arma_order[1] != 1 ||
      !arma_order[2] %in% c(0, 1)) {
    stop(
      "`arma_order` must currently be either c(1, 1) or c(1, 0)."
    )
  }

  if (any(dat$t < 1) || any(dat$t != as.integer(dat$t))) {
    stop("`t` must consist of positive integer time indices.")
  }

  # ==========================================================================

  # 2. Initial response, residuals, and daily aggregation

  # ==========================================================================

  y <- dat[[resp_var]]

  dat$resid <- as.numeric(stats::residuals(mod))
  dat$resid_initial <- dat$resid

  make_daily_summary <- function(data) {


    out <- data |>
      dplyr::group_by(t) |>
      dplyr::summarise(
        w = mean(resid),
        S2 = stats::var(resid),
        n = dplyr::n(),
        .groups = "drop"
      )

    out <- dplyr::left_join(
      tibble::tibble(t = seq_len(max(out$t))),
      out,
      by = "t"
    )

    out$n <- tidyr::replace_na(out$n, 0)

    out


  }

  df_w <- make_daily_summary(dat)
  df_w_initial <- df_w

  T_len <- nrow(df_w)
  n_t <- df_w$n
  observed_days <- n_t > 0

  # ==========================================================================

  # 3. Helper: estimate individual-level variance

  # ==========================================================================

  estimate_var_eps <- function(df_w) {


    df_w_veps_est <- df_w |>
      dplyr::filter(n > 1)

    denom <- sum(df_w_veps_est$n) - nrow(df_w_veps_est)

    if (nrow(df_w_veps_est) == 0 || denom <= 0) {
      stop(
        "Unable to estimate the individual-level variance. ",
        "At least one time point must contain more than one observation."
      )
    }

    sum(
      (df_w_veps_est$n - 1) * df_w_veps_est$S2,
      na.rm = TRUE
    ) / denom


  }

  # ==========================================================================

  # 4. Helper: fit the aggregate ARMA process

  # ==========================================================================

  fit_aggregate_arma <- function(w) {


    fit <- stats::arima(
      w,
      order = c(arma_order[1], 0, arma_order[2]),
      include.mean = FALSE,
      method = "ML"
    )

    phi <- if ("ar1" %in% names(fit$coef)) {
      unname(fit$coef["ar1"])
    } else {
      0
    }

    theta_b <- if ("ma1" %in% names(fit$coef)) {
      unname(fit$coef["ma1"])
    } else {
      0
    }

    list(
      fit = fit,
      phi = phi,
      theta_b = theta_b,
      var_b = fit$sigma2
    )


  }

  # ==========================================================================

  # 5. Initial ARMA estimates

  # ==========================================================================

  arma_est <- fit_aggregate_arma(df_w$w)

  phi_hat <- arma_est$phi
  theta_b_hat <- arma_est$theta_b
  var_b_hat <- arma_est$var_b

  # ==========================================================================

  # 6. Design matrices

  # ==========================================================================

  X_f <- stats::model.matrix(mod)
  X_f_m <- as.matrix(X_f)

  X_t <- f_Xt_sparse(n_t)

  # ==========================================================================

  # 7. Model-specific initialization

  # ==========================================================================

  if (is_mixed) {


    X_r <- lme4::getME(mod, "Z")

    X_r_sparse <- Matrix::Matrix(
      X_r,
      sparse = TRUE
    )

    var_re <- stats::sigma(mod)^2 *
      Matrix::tcrossprod(
        lme4::getME(mod, "Lambda")
      )


  } else {


    X_r <- NULL
    X_r_sparse <- NULL
    var_re <- NULL


  }

  l1 <- as.numeric(stats::logLik(mod))

  # ==========================================================================

  # 8. Initial individual-level variance and Wong transformation

  # ==========================================================================

  var_eps_hat <- estimate_var_eps(df_w)

  var_a_hat <- var_eps_hat *
    mean((n_t[observed_days])^(-1))

  wong_est <- wong(
    ar_hat = phi_hat,
    ma_b_hat = theta_b_hat,
    v_b_hat = var_b_hat,
    v_a_hat = var_a_hat,
    theta_free = arma_order[2] == 1
  )

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]

  # ==========================================================================

  # 9. Initial latent time-series covariance and estimate

  # ==========================================================================

  G_hat <- fG(
    size = T_len,
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

  u_hat <- as.numeric(u_hat)
  u_hat_initial <- u_hat

  # ==========================================================================

  # 10. Helper: calculate inverse covariance components

  # ==========================================================================

  calculate_ts_inverse <- function(G_hat, X_t, var_eps_hat) {


    G_inv <- chol2inv(chol(G_hat))

    Gplus_inv <- chol2inv(
      chol(
        G_inv +
          Matrix::crossprod(X_t) / var_eps_hat
      )
    )

    list(
      G_inv = G_inv,
      Gplus_inv = Gplus_inv
    )


  }

  ts_inverse <- calculate_ts_inverse(
    G_hat = G_hat,
    X_t = X_t,
    var_eps_hat = var_eps_hat
  )

  G_inv <- ts_inverse$G_inv
  Gplus_inv <- ts_inverse$Gplus_inv

  # ==========================================================================

  # 11. Initial level-2 likelihood

  # ==========================================================================

  calculate_l2_fixed <- function(resid,
                                 u_hat,
                                 G_inv,
                                 phi_hat,
                                 var_eta_hat,
                                 T_len) {


    as.numeric(
      -T_len / 2 * log(2 * pi) -
        T_len / 2 * log(var_eta_hat) +
        1 / 2 * log(1 - phi_hat^2) -
        1 / 2 * crossprod(u_hat, G_inv %*% u_hat)
    )


  }

  calculate_l2_mixed <- function(u_hat,
                                 G_hat,
                                 observed_days,
                                 T_len) {


    G_obs <- G_hat[
      observed_days,
      observed_days,
      drop = FALSE
    ]

    G_obs_inv <- solve(G_obs)

    as.numeric(
      -1 / 2 * (
        T_len * log(2 * pi) +
          as.numeric(
            determinant(
              G_obs,
              logarithm = TRUE
            )$modulus
          ) +
          crossprod(
            u_hat,
            G_obs_inv %*% u_hat
          )
      )
    )


  }

  if (is_fixed) {


    Oinvresid <-
      dat$resid / var_eps_hat -
      var_eps_hat^(-2) *
      as.vector(
        X_t %*%
          (
            Gplus_inv %*%
              Matrix::crossprod(X_t, dat$resid)
          )
      )

    u_hat <- as.numeric(
      G_hat %*%
        Matrix::crossprod(X_t, Oinvresid)
    )

    l2 <- calculate_l2_fixed(
      resid = dat$resid,
      u_hat = u_hat,
      G_inv = G_inv,
      phi_hat = phi_hat,
      var_eta_hat = var_eta_hat,
      T_len = T_len
    )


  } else {


    l2 <- calculate_l2_mixed(
      u_hat = u_hat,
      G_hat = G_hat,
      observed_days = observed_days,
      T_len = T_len
    )


  }

  l_new <- l1 + l2
  l_old <- l_new + 2 * tol

  # ==========================================================================

  # 12. Initialize tracking

  # ==========================================================================

  l1_track <- l1
  l2_track <- l2
  l_track <- l_new

  phi_track <- phi_hat
  theta_track <- theta_hat
  var_eta_track <- var_eta_hat
  var_eps_track <- var_eps_hat

  updates <- 0

  # ==========================================================================

  # 13. REMARCS iteration

  # ==========================================================================

  while (
    abs(l_old - l_new) > tol &&
    updates < max_updates
  ) {


    updates <- updates + 1
    l_old <- l_new


    # ------------------------------------------------------------------------
    # 13a. Refresh time indexing
    # ------------------------------------------------------------------------

    n_t <- df_w$n
    observed_days <- n_t > 0

    X_t <- f_Xt_sparse(n_t)


    # ------------------------------------------------------------------------
    # 13b. Remove current latent time-series estimate
    # ------------------------------------------------------------------------

    dat$less_u_hat <- y -
      as.numeric(
        X_t[, observed_days, drop = FALSE] %*%
          u_hat
      )


    # ------------------------------------------------------------------------
    # 13c. Refit regression model for level-1 likelihood
    # ------------------------------------------------------------------------

    new_formula <- stats::formula(mod)
    new_formula[[2]] <- quote(less_u_hat)

    if (is_fixed) {

      reg_new <- stats::lm(
        formula = new_formula,
        data = dat
      )

    } else {

      reg_new <- lme4::lmer(
        formula = new_formula,
        data = dat
      )
    }

    l1 <- as.numeric(stats::logLik(reg_new))


    # ------------------------------------------------------------------------
    # 13d. Time-series inverse covariance calculations
    # ------------------------------------------------------------------------

    ts_inverse <- calculate_ts_inverse(
      G_hat = G_hat,
      X_t = X_t,
      var_eps_hat = var_eps_hat
    )

    G_inv <- ts_inverse$G_inv
    Gplus_inv <- ts_inverse$Gplus_inv


    # ------------------------------------------------------------------------
    # 13e. Update regression parameters
    # ------------------------------------------------------------------------

    if (is_fixed) {

      Oinvy <-
        y / var_eps_hat -
        var_eps_hat^(-2) *
        as.vector(
          X_t %*%
            (
              Gplus_inv %*%
                Matrix::crossprod(X_t, y)
            )
        )

      OinvXf <-
        X_f_m / var_eps_hat -
        var_eps_hat^(-2) *
        as.matrix(
          X_t %*%
            (
              Gplus_inv %*%
                Matrix::crossprod(X_t, X_f_m)
            )
        )

      var_bet <- solve(
        Matrix::crossprod(
          X_f_m,
          OinvXf
        )
      )

      beta_hat <- as.numeric(
        var_bet %*%
          Matrix::crossprod(X_f_m, Oinvy)
      )

      dat$resid <- as.numeric(
        y -
          X_f_m %*%
          beta_hat
      )

    } else {

      beta_hat_2 <- estimate_beta_gls(
        y = y,
        X_f = X_f_m,
        X_r = X_r_sparse,
        X_t = X_t,
        Sigma_b = var_re,
        Gplus_inv = Gplus_inv,
        sigma_eps2 = var_eps_hat
      )

      beta_hat <- as.numeric(beta_hat_2$beta_hat)

      b_hat_list <- estimate_b(
        y = y,
        X_f = X_f_m,
        X_r = X_r_sparse,
        X_t = X_t,
        beta_hat = beta_hat,
        Sigma_b = var_re,
        Gplus_inv = Gplus_inv,
        sigma_eps2 = var_eps_hat
      )

      b_hat <- as.numeric(b_hat_list$b_hat)
      B <- b_hat_list$B

      B_inv <- solve(B)

      var_re_update <- update_Sigma_b(
        b_hat = b_hat,
        B_inv = B_inv,
        mod = mod
      )

      var_re <- var_re_update$matrix

      dat$resid <- as.numeric(
        y -
          X_f_m %*% beta_hat -
          X_r_sparse %*% b_hat
      )
    }


    # ------------------------------------------------------------------------
    # 13f. Update daily residual summaries
    # ------------------------------------------------------------------------

    df_w <- make_daily_summary(dat)

    n_t <- df_w$n
    observed_days <- n_t > 0

    X_t <- f_Xt_sparse(n_t)


    # ------------------------------------------------------------------------
    # 13g. Re-estimate aggregate ARMA process
    # ------------------------------------------------------------------------

    arma_est <- fit_aggregate_arma(df_w$w)

    phi_hat <- arma_est$phi
    theta_b_hat <- arma_est$theta_b
    var_b_hat <- arma_est$var_b


    # ------------------------------------------------------------------------
    # 13h. Update individual-level variance
    # ------------------------------------------------------------------------

    var_eps_hat <- estimate_var_eps(df_w)

    var_a_hat <- var_eps_hat *
      mean((n_t[observed_days])^(-1))


    # ------------------------------------------------------------------------
    # 13i. Wong transformation
    # ------------------------------------------------------------------------

    wong_est <- wong(
      ar_hat = phi_hat,
      ma_b_hat = theta_b_hat,
      v_b_hat = var_b_hat,
      v_a_hat = var_a_hat,
      theta_free = arma_order[2] == 1
    )

    theta_hat <- wong_est[1]
    var_eta_hat <- wong_est[2]


    # ------------------------------------------------------------------------
    # 13j. Update latent time-series covariance
    # ------------------------------------------------------------------------

    G_hat <- fG(
      size = T_len,
      eta_variance = var_eta_hat,
      ar = phi_hat,
      ma = theta_hat
    )


    # ------------------------------------------------------------------------
    # 13k. Update inverse covariance matrices
    # ------------------------------------------------------------------------

    ts_inverse <- calculate_ts_inverse(
      G_hat = G_hat,
      X_t = X_t,
      var_eps_hat = var_eps_hat
    )

    G_inv <- ts_inverse$G_inv
    Gplus_inv <- ts_inverse$Gplus_inv


    # ------------------------------------------------------------------------
    # 13l. Update latent time series
    # ------------------------------------------------------------------------

    if (is_fixed) {

      Oinvresid <-
        dat$resid / var_eps_hat -
        var_eps_hat^(-2) *
        as.vector(
          X_t %*%
            (
              Gplus_inv %*%
                Matrix::crossprod(X_t, dat$resid)
            )
        )

      u_hat <- as.numeric(
        G_hat %*%
          Matrix::crossprod(X_t, Oinvresid)
      )

    } else {

      u_hat <- as.numeric(
        f_u_hat(
          n_t,
          df_w$w,
          G_hat,
          var_eps_hat
        )
      )
    }


    # ------------------------------------------------------------------------
    # 13m. Update level-2 likelihood
    # ------------------------------------------------------------------------

    if (is_fixed) {

      l2 <- calculate_l2_fixed(
        resid = dat$resid,
        u_hat = u_hat,
        G_inv = G_inv,
        phi_hat = phi_hat,
        var_eta_hat = var_eta_hat,
        T_len = T_len
      )

    } else {

      l2 <- calculate_l2_mixed(
        u_hat = u_hat,
        G_hat = G_hat,
        observed_days = observed_days,
        T_len = T_len
      )
    }

    l_new <- l1 + l2


    # ------------------------------------------------------------------------
    # 13n. Track estimates
    # ------------------------------------------------------------------------

    l1_track <- append(l1_track, l1)
    l2_track <- append(l2_track, l2)
    l_track <- append(l_track, l_new)

    phi_track <- append(phi_track, phi_hat)
    theta_track <- append(theta_track, theta_hat)
    var_eta_track <- append(var_eta_track, var_eta_hat)
    var_eps_track <- append(var_eps_track, var_eps_hat)


  }

  # ==========================================================================

  # 14. Final time-series parameter covariance

  # ==========================================================================

  S_o <- Sigma_o_varyingn(
    phi_hat,
    theta_b_hat,
    var_eps_hat,
    var_b_hat,
    n_t
  )

  wong_est <- wong(
    ar_hat = phi_hat,
    ma_b_hat = theta_b_hat,
    v_b_hat = var_b_hat,
    v_a_hat = var_a_hat,
    theta_free = arma_order[2] == 1
  )

  theta_hat <- wong_est[1]
  var_eta_hat <- wong_est[2]

  if (arma_order[2] == 1) {


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

    S_f <- GGG %*%
      S_o %*%
      Matrix::t(GGG)

    S_l <- S_f[
      c(1, 5, 6),
      c(1, 5, 6),
      drop = FALSE
    ]

    S_l_plot <- as.matrix(S_l / T_len)

    ts.pars <- c(
      phi = phi_hat,
      theta = theta_hat
    )

    ts.se <- sqrt(diag(S_l_plot)[1:2])


  } else {


    GG <- matrix(
      0,
      nrow = 3,
      ncol = 4
    )

    GG[1, 1] <- 1
    GG[2, 4] <- 1

    GG[3, 1] <- -2 * phi_hat * var_a_hat
    GG[3, 2] <- 2 * theta_b_hat * var_b_hat
    GG[3, 3] <- 1 + theta_b_hat^2
    GG[3, 4] <- -(1 + phi_hat^2)

    S_f <- GG %*%
      S_o %*%
      t(GG)

    S_l_plot <- as.matrix(S_f / T_len)

    ts.pars <- c(phi = phi_hat)
    ts.se <- sqrt(diag(S_l_plot)[1])


  }

  # ==========================================================================

  # 15. Time-series parameter p-values

  # ==========================================================================

  ts.pvals <- 2 * (
    1 -
      stats::pnorm(
        abs(ts.pars / ts.se)
      )
  )

  ts.pars <- cbind(
    Estimate = ts.pars,
    Std_Error = ts.se,
    P_Value = ts.pvals
  )

  # ==========================================================================

  # 16. Final regression covariance

  # ==========================================================================

  if (is_fixed) {


    # Recalculate with final covariance estimates.

    OinvXf <-
      X_f_m / var_eps_hat -
      var_eps_hat^(-2) *
      as.matrix(
        X_t %*%
          (
            Gplus_inv %*%
              Matrix::crossprod(X_t, X_f_m)
          )
      )

    var_bet <- solve(
      Matrix::crossprod(
        X_f_m,
        OinvXf
      )
    )


  } else {


    var_bet <- as.matrix(
      var_beta(
        X_f = X_f_m,
        X_r = X_r_sparse,
        X_t = X_t,
        Sigma_b = var_re,
        Gplus_inv = Gplus_inv,
        sigma_eps2 = var_eps_hat
      )
    )


  }

  # ==========================================================================

  # 17. Final regression estimates and p-values

  # ==========================================================================

  beta_se <- sqrt(diag(var_bet))

  beta_pvals <- 2 * (
    1 -
      stats::pnorm(
        abs(beta_hat / beta_se)
      )
  )

  beta_out <- cbind(
    Estimate = beta_hat,
    Std_Error = beta_se,
    P_Value = beta_pvals
  )

  rownames(beta_out) <- colnames(X_f_m)

  # ==========================================================================

  # 18. Final fitted values

  # ==========================================================================

  reg_fitted <- as.numeric(
    X_f_m %*%
      beta_hat
  )

  ts_fitted <- as.numeric(
    X_t[, observed_days, drop = FALSE] %*%
      u_hat
  )

  fitted.mat <- data.frame(
    Day = dat$t,
    Reg = reg_fitted,
    TS = ts_fitted
  )

  # ==========================================================================

  # 19. Final residual diagnostic

  # ==========================================================================

  dat$final_resid <-
    y -
    fitted.mat$Reg -
    fitted.mat$TS

  final_resid_daily <- dat |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      mean_resid = mean(final_resid),
      .groups = "drop"
    )

  final_resid_daily <- dplyr::left_join(
    tibble::tibble(
      t = seq_len(max(final_resid_daily$t))
    ),
    final_resid_daily,
    by = "t"
  )

  # ==========================================================================

  # 20. Count estimated parameters for AIC

  # ==========================================================================

  n_beta_pars <- length(beta_hat)

  # phi + var_eta + var_eps, plus theta for ARMA(1,1)

  n_ts_pars <- if (arma_order[2] == 1) 4 else 3

  n_re_pars <- if (is_mixed) {
    length(var_re_update$variances)
  } else {
    0
  }

  n_parameters <- n_beta_pars +
    n_ts_pars +
    n_re_pars

  AIC <- 2 * n_parameters -
    2 * l_new

  # ==========================================================================

  # 21. Build fixed-effects output

  # ==========================================================================

  if (is_fixed) {


    res.list <- list(
      Beta_Pars = beta_out,
      TS_Pars = ts.pars,
      `log(likelihood)` = l_new,
      AIC = AIC,
      Fitted_Values = fitted.mat,
      `Individual Variance` = var_eps_hat,
      `TS Variance` = var_eta_hat,
      `Fixed Covariance` = var_bet,
      `TS Covariance` = S_l_plot,
      Iterations = updates,
      `Likelihood Track` = l_track
    )

    return(res.list)


  }

  # ==========================================================================

  # 22. Generic reconstruction of random effects

  # ==========================================================================

  ranef_template <- lme4::ranef(mod)

  ranef_final <- vector(
    mode = "list",
    length = length(ranef_template)
  )

  names(ranef_final) <- names(ranef_template)

  start <- 1

  for (grp in names(ranef_template)) {


    template <- ranef_template[[grp]]

    n_values <- nrow(template) *
      ncol(template)

    end <- start + n_values - 1

    if (end > length(b_hat)) {
      stop(
        "The length of `b_hat` is incompatible with the random-effects ",
        "structure of `mod`."
      )
    }

    ranef_final[[grp]] <- matrix(
      b_hat[start:end],
      nrow = nrow(template),
      ncol = ncol(template),
      dimnames = dimnames(template)
    )

    start <- end + 1


  }

  if (start - 1 != length(b_hat)) {
    warning(
      "Not all elements of `b_hat` were used when reconstructing ",
      "the random effects."
    )
  }

  class(ranef_final) <- class(ranef_template)

  # ==========================================================================

  # 23. Build mixed-effects output

  # ==========================================================================

  res.list <- list(
    Beta_Pars = beta_out,
    `Random Effects` = ranef_final,
    `Random Effect Variances` = var_re_update$variances,
    TS_Pars = ts.pars,
    `level 1 log(likelihood)` = l1_track,
    `level 2 log(likelihood)` = l2_track,
    `Likelihood Track` = l_track,
    AIC = AIC,
    Fitted_Values = fitted.mat,
    `Time Series Variance` = var_eta_hat,
    `Individual Variance` = var_eps_hat,
    `Initial Time Series Estimate` = u_hat_initial,
    Iterations = updates,
    `TS Covariance` = S_l_plot,
    `Fixed Covariance` = var_bet
  )

  return(res.list)
}
