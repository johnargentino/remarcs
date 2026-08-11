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
RCSfixed.fit <- function(mod, dat, resp_var, tol = 0.1)
{
  resid = residuals(mod)
  dat$resid = resid

  df_w = dat |>
    group_by(t) |>
    summarize(w = mean(resid), S2 = var(resid), n = n())

  plot.orig = ggplot(data = dat) +
    geom_point(mapping = aes(x = t, y = dat[,resp_var]))

  df_w = left_join(tibble(t = 1:max(df_w$t)), df_w)
  acf(df_w$w, na.action = na.pass, main = "ACF of average residuals")

  #Generate initial ARMA(1,1) estimates
  arma_w = arima(df_w$w, order = c(1,0,1), include.mean = FALSE)
  phi_hat = arma_w$coef["ar1"]
  theta_b_hat = arma_w$coef["ma1"]
  var_b_hat = arma_w$sigma2


  df_w$n = replace_na(df_w$n, 0)
  df_w_initial = df_w

  X_f = model.matrix(mod)

  X_f_m = X_f

  l1 = logLik(mod) |> as.numeric()
  dat$resid_initial = dat$resid

  #Get initial estimates
  n_t = df_w$n
  X_t = f_Xt(n_t)
  V = f_V(n_t)
  T = length(n_t)
  N = diag(n_t)

  df_w_veps_est = df_w |>
    filter(n > 1)

  var_eps_hat = sum((df_w_veps_est$n - 1) * df_w_veps_est$S2) / (sum(df_w_veps_est$n) - length(df_w_veps_est$n))

  var_a_hat = var_eps_hat * mean((n_t[n_t != 0])^{-1})
  wong_est = wong(ar_hat=phi_hat, ma_b_hat=theta_b_hat,v_b_hat=var_b_hat,v_a_hat=var_a_hat)
  theta_hat = wong_est[1]
  var_eta_hat = wong_est[2]
  #cat("wong_est = ", wong_est, "\n")
  G_hat = fG(size = T, eta_variance = var_eta_hat, ar = phi_hat, ma = theta_hat)

  u_hat = f_u_hat(n_t, df_w$w, G_hat, var_eps_hat)

  l2 = -1/2 * (T * log(2 * pi) + (determinant(G_hat[n_t != 0, n_t != 0], logarithm = TRUE)$modulus |> as.numeric()) + as.vector(t(u_hat) %*% solve(G_hat[n_t != 0, n_t != 0]) %*% u_hat))
  #cat("l2 exp = ", -1/2*(as.vector(t(u_hat) %*% solve(G_hat[n_t != 0, n_t != 0]) %*% u_hat)), "\n")
  #cat("l2 const = ", -1/2 * (T * log(2 * pi) + (determinant(G_hat[n_t != 0, n_t != 0], logarithm = TRUE)$modulus |> as.numeric())), "\n")
  #cat("l2 = ", l2, "\n")
  l_new = l2 + l1 #Current likelihood
  l_old = l_new + 2 * tol #Why +3?
  l1_track = l1
  l2_track = l2
  l_track = l_old
  phi_track = phi_hat
  theta_track = theta_hat
  var_eta_track = var_eta_hat
  var_eps_track = var_eps_hat


  while (abs(l_old-l_new) > tol){
    l_old = l_new
    NGN_eigen = eigen(N^(1/2) %*% G_hat %*% N ^ (1/2))
    NGN_evec = NGN_eigen$vectors
    NGN_eval = NGN_eigen$values
    L = diag(NGN_eval)
    V = f_V(n_t)
    dat$less_u_hat = dat[[resp_var]] - as.vector(X_t[,n_t != 0] %*% u_hat)


    lm_new = update(mod,less_u_hat~.,data=dat)

    l1 = lm_new |> logLik() |> as.numeric() #This is the likelihood of the fitted model after subtracting the time series
    cat("new model log(lik) = ", l1, "\n")


    beta_hat1 = t(X_f_m) %*% X_t[,n_t != 0] %*% diag(sqrt(1/(n_t[n_t != 0]))) %*% NGN_evec[n_t != 0, n_t != 0] %*% diag(sqrt(1 / (NGN_eval + var_eps_hat)[n_t != 0]))
    beta_hat1 = beta_hat1 %*% t(beta_hat1)

    beta_hat2 = t(X_f_m) %*% V
    tbeta_hat2 = beta_hat2 |> as.matrix() |> t() |> Matrix()
    beta_hat2 = beta_hat2 %*% tbeta_hat2 / var_eps_hat

    beta_hat = solve(beta_hat1 + beta_hat2)
    beta_hat3 = t(X_f_m) %*% X_t[,n_t != 0] %*% diag(sqrt(1/(n_t[n_t != 0]))) %*% NGN_evec[n_t != 0,n_t != 0] %*% diag(1 / (NGN_eval + var_eps_hat)[n_t != 0]) %*% t(NGN_evec[n_t != 0,n_t != 0]) %*% diag(sqrt(1/(n_t)[n_t != 0])) %*% t(X_t[,n_t != 0]) %*% dat[[resp_var]]

    beta_hat4 = Matrix::t(V) %*% dat[[resp_var]]

    beta_hat4 = t(X_f_m) %*% V %*% beta_hat4 / var_eps_hat

    beta_hat = beta_hat %*% (beta_hat3 + beta_hat4)


    dat$resid = as.vector(dat[[resp_var]]- X_f_m %*% beta_hat)

    df_w = dat |>
      group_by(t) |>
      summarize(w = mean(resid), S2 = var(resid), n = n())

    df_w = left_join(tibble(t = 1:max(df_w$t)), df_w)
    df_w$n = replace_na(df_w$n, 0)

    arma_w = arima(df_w$w, order = c(1,0,1), include.mean = FALSE)
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
    G_hat = fG(size = T, eta_variance = var_eta_hat, ar = phi_hat, ma = theta_hat)
    u_hat = f_u_hat(n_t, df_w$w, G_hat, var_eps_hat)
    l2 = -1/2 * (T * log(2 * pi) + (determinant(G_hat[n_t != 0, n_t != 0], logarithm = TRUE)$modulus |> as.numeric()) + as.vector(t(u_hat) %*% solve(G_hat[n_t != 0, n_t != 0]) %*% u_hat))
    cat("l2 = ", l2, "\n")
    l_new = l2 + l1
    l1_track = l1_track |> append(l1)
    l2_track = l2_track |> append(l2)
    l_track = l_track |> append(l_new) #overall likelihood
    cat("likelihood sequence for regression model: \n")

    cat("likelihood sequence for time series model: \n")

    phi_track = phi_track |> append(phi_hat)
    theta_track = theta_track |> append(theta_hat)
    var_eta_track = var_eta_track |> append(var_eta_hat)
    var_eps_track = var_eps_track |> append(var_eps_hat)
  }

  S_o = Sigma_o_varyingn(phi_hat, theta_b_hat, var_eps_hat, var_b_hat, n_t)

  wong_est = wong(phi_hat, theta_b_hat, var_b_hat, var_eps_hat/mean(n_t[n_t != 0]))
  theta_hat = wong_est[1]
  var_eta_hat = wong_est[2]
  var_eta_est = var_eta_hat

  H1 = dHf_dtau1(theta_hat, var_eta_hat)
  H0 = dHf_dtau0(phi_hat, theta_b_hat, var_a_hat, var_b_hat)
  GGG = -solve(H1) %*% H0
  S_f = GGG %*% S_o %*% t(GGG)
  S_l = S_f[c(1,5,6),c(1,5,6)]
  S_l_plot = S_l / T
  vec_est = c(phi_hat, theta_hat, var_eta_hat)
  vec_true = c(phi_hat, theta_hat, var_eta_hat)
  vec_name = c("phi", "theta", "var_eta")

  #plot_ellipsoid(name_vec = vec_name, est_vec = vec_est, true_vec = vec_true, Sigma = S_l_plot)

  XftXf= (t(X_f_m) %*% X_f_m)
  XftXf_inv = solve(XftXf)

  Sigma_beta3_term1 = t(X_f_m) %*% X_t[, n_t != 0] %*% diag(sqrt(1 / n_t[n_t != 0])) %*% NGN_evec[n_t != 0,n_t != 0] %*% diag(sqrt(1/(NGN_eval[n_t != 0] + var_eps_hat)))

  Sigma_beta3_term1 = Sigma_beta3_term1 %*% Matrix::t(Sigma_beta3_term1)

  Sigma_beta3_term2 = t(X_f_m) %*% V
  Sigma_beta3_term2 = Sigma_beta3_term2 %*% Matrix::t(Sigma_beta3_term2) / var_eps_hat

  Sigma_beta3 = solve(Sigma_beta3_term1 + Sigma_beta3_term2)

  beta_hat = cbind(beta_hat, sqrt(diag(Sigma_beta3)))
  pvals=2*(1-pnorm(abs(beta_hat[,1]/beta_hat[,2])))
  beta_hat=cbind(beta_hat,pvals)

  AIC=2*(nrow(beta_hat)+3)-2*l_new

  #armasim = arima.sim(model = list(ar = phi_hat, ma = theta_hat), sd = sqrt(var_eta_hat), n = T)

  #plot.1=ggplot() +
  #  geom_point(mapping = aes(x = unique(dat$t), y = u_hat), col = "blue")

  #plot.2= ggplot() +
  #  geom_point(mapping = aes(x = 1:T, y = armasim))

  #plot.3=ggplot() +
  #  geom_point(data = dat, mapping = aes(x = t, y = dat[,resp_var]), size = .1) +
  #  geom_point(mapping = aes(x = dat$t, y = X_t[,n_t != 0] %*% u_hat + X_f %*% beta_hat[,1]), col = "blue", size = .1)

  fitted.mat=data.frame(Day=dat$t,Reg=X_f %*% beta_hat[,1],TS=X_t[,n_t != 0] %*% u_hat)



  ts.pars=c(phi_hat,theta_hat,var_eta_hat)
  names(ts.pars)=c("phi","theta","var_eta")
  res.list=list(beta_hat,ts.pars,l_new, AIC,fitted.mat)
  names(res.list)=c("Beta_Pars","TS_Pars","log(likelihood)", "AIC","Fitted_Values")

  return(res.list)
}
