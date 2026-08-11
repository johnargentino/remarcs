# Internal helper functions for the REMARCS package
#
# These functions are used internally by RCSmixed.fit() and are not exported.

# Construct the observation-by-time design matrix for the latent time series.
# Missing days are represented by zero entries in n_vec.
f_Xt <- function(n_vec) {
  X_t <- matrix(
    0,
    nrow = sum(n_vec),
    ncol = length(n_vec)
  )

  trackrow <- 1L

  for (i in seq_along(n_vec)) {
    if (n_vec[i] > 0) {
      X_t[
        trackrow:(trackrow + n_vec[i] - 1L),
        i
      ] <- 1
      trackrow <- trackrow + n_vec[i]
    }
  }

  X_t
}

equations_arma <- function(vars, ph_hat, theta_b_hat, sigma_a2_hat, sigma_b2_hat) {
  th <- vars[1]
  sigma_eps2 <- vars[2]

  eq1 <- (1 + th^2) * sigma_eps2 + (1 + ph_hat^2) * sigma_a2_hat - (1 + theta_b_hat^2) * sigma_b2_hat
  eq2 <- th * sigma_eps2 - ph_hat * sigma_a2_hat - theta_b_hat * sigma_b2_hat

  return(eq1^2 + eq2^2)
}
# Raise base to exponent, treating negative exponents as zero.
power_f <- function(base, exponent) {
  exponent <- pmax(exponent, 0)
  base^exponent
}


# Helper used to construct the ARMA(1,1) covariance matrix.
outer_function_G <- function(x, y, ar_coef, ma_coef, vv) {
  w <- abs(x - y)

  z <- vv *
    (
      ar_coef +
        ma_coef +
        ar_coef * (ar_coef + ma_coef)^2 / (1 - ar_coef^2)
    ) *
    power_f(ar_coef, w - 1)

  as.vector(z)
}


# Construct the covariance matrix of the latent ARMA(1,1) process.
fG <- function(size, eta_variance, ar, ma) {
  dd <- eta_variance * (
    1 + (ar + ma)^2 / (1 - ar^2)
  )

  z <- outer(
    1:size,
    1:size,
    FUN = function(x, y) {
      outer_function_G(
        x,
        y,
        ar_coef = ar,
        ma_coef = ma,
        vv = eta_variance
      )
    }
  )

  z <- z -
    diag(
      outer_function_G(
        1,
        1,
        ar_coef = ar,
        ma_coef = ma,
        vv = eta_variance
      ),
      size
    ) +
    diag(dd, size)

  z
}


# Asymptotic covariance matrix used for the varying-n case.
Sigma_o_varyingn <- function(ar, ma_b, veps, v_b, n_t) {
  block11 <- (1 + ar * ma_b) / (ar + ma_b)^2 *
    cbind(
      c(
        (1 - ar^2) * (1 + ar * ma_b),
        -(1 - ar^2) * (1 - ma_b^2)
      ),
      c(
        -(1 - ar^2) * (1 - ma_b^2),
        (1 - ma_b^2) * (1 + ar * ma_b)
      )
    )

  block22 <- 2 * v_b^2

  n_t <- n_t[n_t != 0]
  n_bar <- mean(n_t)
  m_hat <- mean(1 / n_t)
  S2_n_inverse <- stats::var(1 / n_t)

  avar_sigma_a2 <- S2_n_inverse * veps^2 +
    2 * veps^2 * m_hat^2 / (n_bar - 1)

  block33 <- avar_sigma_a2

  full <- matrix(0, nrow = 4, ncol = 4)
  full[1:2, 1:2] <- block11
  full[3, 3] <- block22
  full[4, 4] <- block33

  full
}


# Recover the latent ARMA(1,1) parameters from the observed noisy process.
wong <- function(
    ar_hat,
    ma_b_hat,
    v_b_hat,
    v_a_hat,
    theta_free = TRUE) {

  theta_b_hat <- ma_b_hat
  ph_hat <- ar_hat
  sigma_a2_hat <- v_a_hat
  sigma_b2_hat <- v_b_hat

  if (!theta_free) {
    theta_hat <- 0
    var_eta_hat <- max(
      (1 + theta_b_hat^2) * sigma_b2_hat -
        (1 + ph_hat^2) * sigma_a2_hat,
      0
    )

    return(c(theta_hat, var_eta_hat))
  }

  # Initial guesses
  aa <- theta_b_hat * sigma_b2_hat +
    ph_hat * sigma_a2_hat

  bb <- (1 + ph_hat^2) * sigma_a2_hat -
    (1 + theta_b_hat^2) * sigma_b2_hat

  cc <- aa

  start_ma <- 1 / (2 * aa) *
    (-bb - sqrt(bb^2 - 4 * aa * cc))

  start_veta <- (
    theta_b_hat * sigma_b2_hat +
      ph_hat * sigma_a2_hat
  ) / start_ma

  if (!is.na(start_ma) && !is.na(start_veta)) {
    return(
      c(
        start_ma,
        start_veta,
        !is.na(start_ma) && !is.na(start_veta)
      )
    )
  }

  start <- c(theta_b_hat, sigma_b2_hat)

  # Constraints
  lower_bounds <- c(-0.9999, 1e-6)
  upper_bounds <- c(0.9999, Inf)

  # Solve the system
  result <- BB::BBoptim(
    par = start,
    fn = equations_arma,
    ph_hat = ph_hat,
    theta_b_hat = theta_b_hat,
    sigma_a2_hat = sigma_a2_hat,
    sigma_b2_hat = sigma_b2_hat,
    lower = lower_bounds,
    upper = upper_bounds,
    quiet = TRUE,
    control = list(trace = 0)
  )

  theta_hat <- result$par[1]
  var_eta_hat <- result$par[2]

  c(
    theta_hat,
    var_eta_hat,
    !is.na(start_ma) && !is.na(start_veta)
  )
}


# Derivative of the transformation from observed to latent ARMA parameters.
dHf_dtau0 <- function(phi, theta_b, sigma_a2, sigma_b2) {
  M <- matrix(0, nrow = 6, ncol = 4)

  # Identity block (negative)
  diag(M[1:4, 1:4]) <- -1

  # Row 5
  M[5, 1] <- 2 * phi * sigma_a2
  M[5, 2] <- -2 * theta_b * sigma_b2
  M[5, 3] <- -(1 + theta_b^2)
  M[5, 4] <- 1 + phi^2

  # Row 6
  M[6, 1] <- -sigma_a2
  M[6, 2] <- -sigma_b2
  M[6, 3] <- -theta_b
  M[6, 4] <- -phi

  M
}


# Derivative of the transformation from latent to observed ARMA parameters.
dHf_dtau1 <- function(theta, sigma_eta2) {
  M <- matrix(0, nrow = 6, ncol = 6)

  # Identity block
  diag(M[1:4, 1:4]) <- 1

  # Row 5
  M[5, 5] <- 2 * theta * sigma_eta2
  M[5, 6] <- 1 + theta^2

  # Row 6
  M[6, 5] <- sigma_eta2
  M[6, 6] <- theta

  M
}


# Construct the within-day contrast matrix.
f_V <- function(n_t) {
  T <- length(n_t)
  V_list <- list()
  vlistcount <- 1L

  for (nt in seq_len(T)) {
    nrow <- n_t[nt]

    if (nrow != 0) {
      c_space <- rep(1, nrow)
      null_basis <- MASS::Null(c_space)
      ortho_basis <-qr.Q(qr(null_basis))

      V_list[[vlistcount]] <- ortho_basis
      vlistcount <- vlistcount + 1L
    }
  }

  V <- Matrix::bdiag(V_list)
  Matrix::Matrix(V, sparse = TRUE)
}


# Estimate the conditional expectation of the latent time-series process.
f_u_hat <- function(n_t, residuals, Gamma, var_epsilon) {
  X_t <- f_Xt(n_t)
  sqN <- diag(sqrt(n_t))

  NGN_eigen <- eigen(sqN %*% Gamma %*% sqN)

  NGN_evec <- NGN_eigen$vectors
  NGN_eval <- NGN_eigen$values

  X_t <- Matrix::Matrix(X_t, sparse = TRUE)

  L <- Matrix::Matrix(
    diag(NGN_eval + var_epsilon),
    sparse = TRUE
  )

  u1 <- Gamma %*%
    sqN %*%
    NGN_evec %*%
    solve(L) %*%
    t(NGN_evec) %*%
    sqN

  u1 <- u1[n_t != 0, n_t != 0]
  u1 <- u1 %*% residuals[n_t != 0]

  as.vector(u1)
}

estimate_b <- function(
    y,
    X_f,
    X_r,
    X_t,
    beta_hat,
    Sigma_b,
    Gplus_inv,
    sigma_eps2
) {

  # ---------------------------------------------------------
  # 1. Residual after removing fixed effects
  # ---------------------------------------------------------

  r <- y - X_f %*% beta_hat


  # ---------------------------------------------------------
  # 2. Useful cross-products
  # ---------------------------------------------------------

  XtXr <- Matrix::crossprod(X_t, X_r)  # X_t' X_r
  Xtr  <- Matrix::crossprod(X_t, r)    # X_t' r

  XrXr <- Matrix::crossprod(X_r)       # X_r' X_r
  Xrr  <- Matrix::crossprod(X_r, r)    # X_r' r


  # ---------------------------------------------------------
  # 3. Calculate v = X_r' A^{-1} r
  #
  # A^{-1} =
  #   sigma_eps^-2 I
  #   - sigma_eps^-4 X_t Gplus_inv X_t'
  # ---------------------------------------------------------

  c <- 1 / sigma_eps2

  v <- c * Xrr -
    c^2 * Matrix::crossprod(
      XtXr,
      Gplus_inv %*% Xtr
    )


  # ---------------------------------------------------------
  # 4. Calculate K = X_r' A^{-1} X_r
  # ---------------------------------------------------------

  K <- c * XrXr -
    c^2 * Matrix::crossprod(
      XtXr,
      Gplus_inv %*% XtXr
    )


  # ---------------------------------------------------------
  # 5. B = Sigma_b^{-1} + K
  # ---------------------------------------------------------

  Sigma_b_inv <- solve(Sigma_b)

  B <- Sigma_b_inv + K


  # ---------------------------------------------------------
  # 6. BLUP / posterior mean of b
  #
  # b_hat = B^{-1} v
  # ---------------------------------------------------------

  b_hat <- solve(B, v)


  return(list(
    b_hat = b_hat,
    B = B,
    v = v
  ))
}


estimate_beta_gls <- function(
    y,
    X_f,
    X_r,
    X_t,
    Sigma_b,
    Gplus_inv,
    sigma_eps2
) {

  # ---------------------------------------------------------
  # 1. Constants
  # ---------------------------------------------------------

  c <- 1 / sigma_eps2


  # ---------------------------------------------------------
  # 2. Cross-products involving X_t
  # ---------------------------------------------------------

  XtXr <- Matrix::crossprod(X_t, X_r)  # X_t' X_r
  XtXf <- Matrix::crossprod(X_t, X_f)  # X_t' X_f
  Xty  <- Matrix::crossprod(X_t, y)    # X_t' y


  # ---------------------------------------------------------
  # 3. K = X_r' A^{-1} X_r
  # ---------------------------------------------------------

  XrXr <- Matrix::crossprod(X_r)

  K <- c * XrXr -
    c^2 * Matrix::crossprod(
      XtXr,
      Gplus_inv %*% XtXr
    )


  # ---------------------------------------------------------
  # 4. B = Sigma_b^{-1} + K
  # ---------------------------------------------------------

  Sigma_b_inv <- solve(Sigma_b)

  B <- Sigma_b_inv + K


  # ---------------------------------------------------------
  # 5. A^{-1} X_f
  # ---------------------------------------------------------

  Ainv_Xf <- c * X_f -
    c^2 * X_t %*% (
      Gplus_inv %*% XtXf
    )


  # ---------------------------------------------------------
  # 6. A^{-1} y
  # ---------------------------------------------------------

  Ainv_y <- c * y -
    c^2 * X_t %*% (
      Gplus_inv %*% Xty
    )


  # ---------------------------------------------------------
  # 7. X_r' A^{-1} X_f
  # ---------------------------------------------------------

  Xr_Ainv_Xf <- Matrix::crossprod(
    X_r,
    Ainv_Xf
  )


  # ---------------------------------------------------------
  # 8. X_r' A^{-1} y
  # ---------------------------------------------------------

  Xr_Ainv_y <- Matrix::crossprod(
    X_r,
    Ainv_y
  )


  # ---------------------------------------------------------
  # 9. Compute:
  #
  # X_f' Sigma_y^{-1} X_f
  #
  # = X_f' A^{-1} X_f
  #   - X_f' A^{-1} X_r
  #     B^{-1}
  #     X_r' A^{-1} X_f
  # ---------------------------------------------------------

  Xf_Ainv_Xf <- Matrix::crossprod(
    X_f,
    Ainv_Xf
  )

  correction_Xf <- Matrix::t(Xr_Ainv_Xf) %*%
    solve(B, Xr_Ainv_Xf)

  Xt_Siginv_Xf <- Xf_Ainv_Xf -
    correction_Xf


  # ---------------------------------------------------------
  # 10. Compute:
  #
  # X_f' Sigma_y^{-1} y
  #
  # = X_f' A^{-1} y
  #   - X_f' A^{-1} X_r
  #     B^{-1}
  #     X_r' A^{-1} y
  # ---------------------------------------------------------

  Xf_Ainv_y <- Matrix::crossprod(
    X_f,
    Ainv_y
  )

  correction_y <- Matrix::t(Xr_Ainv_Xf) %*%
    solve(B, Xr_Ainv_y)

  Xf_Siginv_y <- Xf_Ainv_y -
    correction_y


  # ---------------------------------------------------------
  # 11. GLS estimate
  # ---------------------------------------------------------

  beta_hat <- solve(
    Xt_Siginv_Xf,
    Xf_Siginv_y
  )


  # ---------------------------------------------------------
  # 12. Return useful quantities
  # ---------------------------------------------------------

  return(list(
    beta_hat = beta_hat,
    B = B,
    Xt_Siginv_Xf = Xt_Siginv_Xf,
    Xf_Siginv_y = Xf_Siginv_y
  ))
}

update_Sigma_b <- function(b_hat, B_inv, mod) {

  # Number of levels in each random-effect grouping factor
  m1 <- nlevels(getME(mod, "flist")$adm1)
  m2 <- nlevels(getME(mod, "flist")$adm2)

  # Split random effects
  b1 <- b_hat[1:m1]
  b2 <- b_hat[(m1 + 1):(m1 + m2)]

  # Conditional covariance blocks
  C11 <- B_inv[1:m1, 1:m1]
  C22 <- B_inv[(m1 + 1):(m1 + m2),
               (m1 + 1):(m1 + m2)]

  # EM variance updates
  sigma2_adm1 <- (
    sum(b1^2) + sum(diag(C11))
  ) / m1

  sigma2_adm2 <- (
    sum(b2^2) + sum(diag(C22))
  ) / m2

  return(c(
    sigma2_adm1 = sigma2_adm1,
    sigma2_adm2 = sigma2_adm2
  ))
}

update_Sigma_b <- function(b_hat, B_inv, mod) {

  # Number of levels in each random-effect grouping factor
  m1 <- nlevels(getME(mod, "flist")$adm1)
  m2 <- nlevels(getME(mod, "flist")$'adm1:adm2')

  # Split random effects
  b1 <- b_hat[1:m1]
  b2 <- b_hat[(m1 + 1):(m1 + m2)]

  # Conditional covariance blocks
  C11 <- B_inv[1:m1, 1:m1]
  C22 <- B_inv[(m1 + 1):(m1 + m2),
               (m1 + 1):(m1 + m2)]

  # EM variance updates
  sigma2_adm1 <- (
    sum(b1^2) + sum(diag(C11))
  ) / m1

  sigma2_adm2 <- (
    sum(b2^2) + sum(diag(C22))
  ) / m2
  var_re2 = diag(rep(sigma2_adm2, m2))
  var_re1 = diag(rep(sigma2_adm1, m1))
  return(list(sigma2 = c(
    sigma2_adm1 = sigma2_adm1,
    sigma2_adm2 = sigma2_adm2
  ), var_re = bdiag(var_re2, var_re1)))
}
