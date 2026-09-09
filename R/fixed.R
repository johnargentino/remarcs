#' Fit the REMARCS model with a fixed-effects regression component
#'
#' Implements the REMARCS algorithm using an initial fixed-effects model
#' and a latent ARMA(1, 1) time-series component.
#'
#' @param mod The initial fixed-effects model of the response variable.
#' @param dat The data used for fitting. The variable for time must be named
#'   `t`.
#' @param resp_var Character string giving the name of the response variable.
#' @param tol Convergence tolerance for the iterative REMARCS algorithm.
#'
#' @return A list containing the estimated regression parameters, time-series
#'   parameters, log-likelihood, AIC, and fitted values.
#'
#' @export
RCSfixed.fit <- function(mod, dat, resp_var, tol = 0.1, max_updates = 10)
{
  y = dat[[resp_var]]
  resid = residuals(mod)
  dat$resid = resid

  df_w = dat |>
    group_by(t) |>
    summarize(w = mean(resid), S2 = var(resid), n = n())



  df_w = left_join(tibble(t = 1:max(df_w$t)), df_w, by = "t")
  # acf(df_w$w, na.action = na.pass, main = "ACF of average residuals")

  #Generate initial ARMA(1,1) estimates
  arma_w = arima(df_w$w, order = c(1,0,1), include.mean = FALSE, method = "ML")
  phi_hat = arma_w$coef["ar1"]
  theta_b_hat = arma_w$coef["ma1"]
  var_b_hat = arma_w$sigma2


  df_w$n = replace_na(df_w$n, 0)
  df_w_initial = df_w
  T_len = nrow(df_w)

  X_f = model.matrix(mod)

  X_f_m = X_f

  l1 = logLik(mod) |> as.numeric()
  dat$resid_initial = dat$resid

  #Get initial estimates
  n_t = df_w$n
  X_t = f_Xt_sparse(n_t)
  V = f_V(n_t)
  # T = length(n_t)
  N = diag(n_t)
  df_w_veps_est = df_w |>
    filter(n > 1)
  var_eps_hat = sum((df_w_veps_est$n - 1) * df_w_veps_est$S2) / (sum(df_w_veps_est$n) - length(df_w_veps_est$n))
  var_a_hat = var_eps_hat * mean((n_t[n_t != 0])^{-1})
  wong_est = wong(ar_hat=phi_hat, ma_b_hat=theta_b_hat,v_b_hat=var_b_hat,v_a_hat=var_a_hat)
  theta_hat = wong_est[1]
  var_eta_hat = wong_est[2]
  G_hat = fG(size = T_len, eta_variance = var_eta_hat, ar = phi_hat, ma = theta_hat)
  G_inv = G_hat |>
    chol() |>
    chol2inv()
  Gplus_inv = G_inv + Matrix::t(X_t) %*% X_t / var_eps_hat
  Gplus_inv = Gplus_inv |>
    chol() |>
    chol2inv()
  Oinvresid = (dat$resid / var_eps_hat) - var_eps_hat ^ (-2) * (X_t %*% (Gplus_inv %*% (Matrix::t(X_t) %*% dat$resid)))
  u_hat = (G_hat %*% (Matrix::t(X_t) %*% Oinvresid)) %>% as.vector()
  l2 = -T_len/2 * log(2 * pi) - T_len / 2 * log(var_eta_hat) + 1 / 2 * log(1 - phi_hat ^ 2) - 1 / 2 * t(u_hat) %*% G_inv %*% u_hat
  l_new = l2 + l1
  l_old = l_new + 2 * tol
  l1_track = l1
  l2_track = l2
  l_track = l_old
  phi_track = phi_hat
  theta_track = theta_hat
  var_eta_track = var_eta_hat
  var_eps_track = var_eps_hat
updates = 1
  while (abs(l_old-l_new) > tol & updates <= max_updates){
    updates = updates + 1
    l_old = l_new
    dat$less_u_hat = dat[[resp_var]] - as.vector(X_t[,n_t != 0] %*% u_hat)
    lm_new = update(mod,less_u_hat~.,data=dat)
    l1 = lm_new |> logLik() |> as.numeric() #This is the likelihood of the fitted model after subtracting the time series
    Oinvy = (dat[[resp_var]] / var_eps_hat) - var_eps_hat ^ (-2) * (X_t %*% (Gplus_inv %*% (Matrix::t(X_t) %*% dat[[resp_var]])))
    OinvXf = (X_f_m / var_eps_hat) - var_eps_hat ^ (-2) * (X_t %*% (Gplus_inv %*% (Matrix::t(X_t) %*% X_f_m)))
    beta_hat = solve(Matrix::t(X_f_m) %*% OinvXf) %*% Matrix::t(X_f_m) %*% Oinvy
    beta_hat = as.numeric(beta_hat)
    dat$resid = as.vector(dat[[resp_var]]- X_f_m %*% beta_hat)

    df_w = dat |>
      group_by(t) |>
      summarize(w = mean(resid), S2 = var(resid), n = n())

    df_w = left_join(tibble(t = 1:max(df_w$t)), df_w, by = "t")
    df_w$n = replace_na(df_w$n, 0)

    arma_w = arima(df_w$w, order = c(1,0,1), include.mean = FALSE, method = "ML")
    phi_hat = arma_w$coef["ar1"]
    theta_b_hat = arma_w$coef["ma1"]
    var_b_hat = arma_w$sigma2


    df_w_veps_est = df_w |>
      filter(n > 1)

    var_eps_hat = sum((df_w_veps_est$n - 1) * df_w_veps_est$S2) / (sum(df_w_veps_est$n) - length(df_w_veps_est$n))
    var_a_hat = var_eps_hat * mean((n_t[n_t != 0])^{-1})

    wong_est = wong(phi_hat, theta_b_hat, var_b_hat, var_a_hat)
    theta_hat = wong_est[1]
    var_eta_hat = wong_est[2]
    G_hat = fG(size = T_len, eta_variance = var_eta_hat, ar = phi_hat, ma = theta_hat)

    G_inv = G_hat |>
      chol() |>
      chol2inv()
    Gplus_inv = G_inv + Matrix::t(X_t) %*% X_t / var_eps_hat
    Gplus_inv = Gplus_inv |>
      chol() |>
      chol2inv()
    Oinvresid = (dat$resid / var_eps_hat) - var_eps_hat ^ (-2) * (X_t %*% (Gplus_inv %*% (Matrix::t(X_t) %*% dat$resid)))
    u_hat = (G_hat %*% (Matrix::t(X_t) %*% Oinvresid)) %>% as.vector()

    l2 = -T_len/2 * log(2 * pi) - T_len / 2 * log(var_eta_hat) + 1 / 2 * log(1 - phi_hat ^ 2) - 1 / 2 * t(u_hat) %*% G_inv %*% u_hat


    l_new = l2 + l1
    l1_track = l1_track |> append(l1)
    l2_track = l2_track |> append(l2)
    l_track = l_track |> append(l_new) #overall likelihood
    # cat("likelihood sequence for regression model: \n")
    #
    # cat("likelihood sequence for time series model: \n")

    phi_track = phi_track |> append(phi_hat)
    theta_track = theta_track |> append(theta_hat)
    var_eta_track = var_eta_track |> append(var_eta_hat)
    var_eps_track = var_eps_track |> append(var_eps_hat)
  }

  S_o = Sigma_o_varyingn(phi_hat, theta_b_hat, var_eps_hat, var_b_hat, n_t)

print(updates)

  H1 = dHf_dtau1(theta_hat, var_eta_hat)
  H0 = dHf_dtau0(phi_hat, theta_b_hat, var_a_hat, var_b_hat)
  GGG = -solve(H1) %*% H0
  S_f = GGG %*% S_o %*% t(GGG)
  S_l = S_f[c(1,5,6),c(1,5,6)]
  S_l_plot = S_l / T_len
  S_l_plot = as.matrix(S_l_plot)
  vec_est = c(phi_hat, theta_hat, var_eta_hat)
  vec_true = c(phi_hat, theta_hat, var_eta_hat)
  vec_name = c("phi", "theta", "var_eta")
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
  OinvXf = (X_f_m / var_eps_hat) - var_eps_hat ^ (-2) * (X_t %*% (Gplus_inv %*% (Matrix::t(X_t) %*% X_f_m)))
  var_bet = solve(Matrix::t(X_f_m) %*% OinvXf)
  var_bet = as.matrix(var_bet)
  beta_hat = cbind(beta_hat, sqrt(diag(var_bet)))
  pvals=2*(1-pnorm(abs(beta_hat[,1]/beta_hat[,2])))
  beta_hat=cbind(beta_hat,pvals)
  # --- REPLACE THE BOTTOM RETURN BLOCK WITH THIS ---
  AIC = 2 * (nrow(beta_hat) + 3) - 2 * l_new

  # FIX: Use the dense matrix X_f_m instead of the tibble X_f, and drop S4 Matrix tracking
  reg_fitted = as.numeric(X_f_m %*% beta_hat[, 1])
  ts_fitted  = as.numeric(X_t %*% u_hat)

  fitted.mat = data.frame(
    Day = dat$t,
    Reg = reg_fitted,
    TS  = ts_fitted
  )

  res.list = list(
    beta_hat,
    ts.pars,
    l_new,
    AIC,
    fitted.mat,
    var_eps_hat,
    var_eta_hat,
    var_bet,
    S_l_plot
  )

  names(res.list) = c("Beta_Pars", "TS_Pars", "log(likelihood)", "AIC", "Fitted_Values",
                      "Individual Variance", "TS Variance", "Fixed Covariance", "TS Covariance")
  return(res.list)

}
