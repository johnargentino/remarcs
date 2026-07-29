#' remarcs: REMARCS Estimation for Repeated Cross-Sectional Data
#'
#' Tools for estimating repeated cross-sectional regression models with
#' latent temporal dependence using the REMARCS algorithm.
#'
#' @importFrom stats acf arima formula logLik model.matrix
#' @importFrom stats na.pass pnorm residuals sigma update var
#'
#' @importFrom dplyr group_by summarize n pull left_join filter
#' @importFrom ggplot2 ggplot geom_point aes
#' @importFrom tidyr replace_na
#' @importFrom tibble tibble
#' @importFrom lme4 getME VarCorr
#'
#' @importFrom Matrix Matrix
#'
#' @keywords internal
"_PACKAGE"