#' Compute block-diagonal smoothing matrix 
#'
#' @param gm value returned by mgcv::gam with fit = FALSE
#'
#' @return list of S (block diagonal smoothing matrix), Scols (number of columns for each smoothing matrix), 
#' bdim (k argument used for each smooth), Snames (names of smoothing parameters). Returns NULL if no smooths in
#' model. 
#' 
#' @description Based on code written by David L Miller in CTMCdive R Package (https://github.com/r-glennie/CTMCdive). 
#' 
#' @export
#' @importFrom Matrix bdiag
make_smoothing_matrix <- function(gm) {
  if (length(gm$smooth) < 0.5) return(NULL)
  S <- list() 
  Sncols <- list()
  bdim <- list()
  smnames <- list()
  for (i in seq_along(gm$smooth)) {
    sm <- gm$smooth[[i]]
    for (j in seq_along(sm$S)) {
      S <- c(S, as(sm$S[[j]], "sparseMatrix"))
      Sncols <- c(Sncols, ncol(sm$S[[j]]))
      bdim <- c(bdim, sm$bs.dim)
      smnames <- c(smnames, attr(sm$sp, "names"))
    }
  }
  Smat <- bdiag(S)
  res <- list(S = Smat, Scols = Sncols, bdim = bdim, Snames = smnames)
  return(res)
}

#' Compute design and smoothing matrices from formulae 
#'
#' @param forms list of formulae for psi and p 
#' @param site_data see fit_occu function 
#' @param visit_data see fit_occu function 
#'
#' @return list of X = fixed effects design matrix, Z = random effects design matrix, 
#'   and S = smoothing matrix for both psi and p (e.g. X_psi)
#' @export
#'
#' @examples
#' @importFrom mgcv gam 
#' @importFrom mgcv predict.gam 
#' @importFrom methods as 
make_matrices <- function(forms, visit_data, site_data) {
  ## occupancy model 
  site_data$psi <- 1:nrow(site_data)
  gam_psi <- gam(forms[["psi"]], data = site_data, method = "REML")
  site_data <- site_data[, psi := NULL]
  # get design matrix 
  Xfull <- predict(gam_psi, newdata = site_data, type = "lpmatrix")
  X_psi <- Xfull[, 1:gam_psi$nsdf, drop = FALSE]
  Z_psi <- Xfull[, -(1:gam_psi$nsdf), drop = FALSE]
  # get smoothing matrix 
  S_psi <- make_smoothing_matrix(gam_psi)
  
  ## detection model 
  visit_data$p <- 1:nrow(visit_data)
  gam_p <- gam(forms[["p"]], data = visit_data, method = "REML")
  visit_data <- visit_data[, p := NULL]
  # get design matrix
  Xfull <- predict(gam_p, newdata = visit_data, type = "lpmatrix")
  X_p <- Xfull[, 1:gam_p$nsdf, drop = FALSE]
  Z_p <- Xfull[, -(1:gam_p$nsdf), drop = FALSE]
  # get smoothing matrix
  S_p <- make_smoothing_matrix(gam_p)
  
  ## results 
  res <- list(X_psi = X_psi, 
              Z_psi = Z_psi, 
              S_psi = S_psi, 
              X_p = X_p, 
              Z_p = Z_p,
              S_p = S_p)
  
  return(res)
}

#' Fit occupancy model 
#'
#' @param forms list of formulae for psi (occupancy) and p (detection) 
#' @param visit_data data.table with row for each visit with a y column for detection record 
#' @param site_data data.table with row for each site 
#' @return fitted occupied model object 
#' @export
#' @importFrom data.table uniqueN
#' @importFrom TMB MakeADFun sdreport normalize
#' @useDynLib occu_tmb
fit_occu <- function(forms, visit_data, site_data, print = TRUE) {
  
  ## DATA 
  # order data 
  site_data <- site_data[order(site_data$site, site_data$occasion),]
  visit_data <- visit_data[order(visit_data$site, visit_data$occasion, visit_data$visit),]
  # number of sites
  nsites <- uniqueN(site_data$site)
  # number of occasions 
  nocc <- uniqueN(site_data$occasion)
  # index of first visit each site x occasion
  ind <- matrix(as.numeric(which(visit_data$visit == 1)), nr = nsites, nc = nocc) - 1
  # total number of detections per site x occasion
  totsite <- matrix(visit_data[, .(totsite = sum(y)), .(site, occasion)]$totsite, nr = nsites, nc = nocc)
  # number of visits per site x occasion
  nvisit <- matrix(visit_data[, .(nvisit = .N), .(site, occasion)]$nvisit, nr = nsites, nc = nocc)
  
  ## MODEL MATRICES
  # name formulae
  names(forms) <- c(as.character(forms[[1]][[2]]),
                    as.character(forms[[2]][[2]]))
  # get model matrices and smoothing matrices 
  mats <- make_matrices(forms, visit_data, site_data)
  
  ## SETUP DATA FOR TMB
  tmb_dat <- list(flag = 1L, 
                  nsites = nsites, 
                  nocc = nocc, 
                  y = visit_data$y, 
                  ind = ind, 
                  totsite = totsite, 
                  nvisit = nvisit, 
                  X_psi = mats$X_psi, 
                  Z_psi = mats$Z_psi, 
                  S_psi = mats$S_psi$S, 
                  S_psi_n = as.integer(mats$S_psi$Scols), 
                  X_p = mats$X_p, 
                  Z_p = mats$Z_p, 
                  S_p = mats$S_p$S, 
                  S_p_n = as.integer(mats$S_p$Scols))
  
  ## PARAMETERS 
  map <- list()
  random <- NULL
  # fixed effects 
  beta_psi <- rep(0, ncol(tmb_dat$X_psi))
  beta_p <- rep(0, ncol(tmb_dat$X_p))
  # random effects 
  # occupancy
  if (is.null(mats$S_psi)) {
    z_psi <- 0 
    log_lambda_psi <- 0
    tmb_dat$S_psi <- as(matrix(0, 1, 1), "sparseMatrix") 
    tmb_dat$S_psi_n <- -1
    map <- c(map, list(z_psi = as.factor(NA), log_lambda_psi = as.factor(NA)))
  } else {
    z_psi <- rep(0, ncol(tmb_dat$Z_psi)) 
    log_lambda_psi <- rep(0, length(tmb_dat$S_psi_n)) 
    random <- c(random, "z_psi")
  }
  # detection 
  if (is.null(mats$S_p)) {
    z_p <- 0 
    log_lambda_p <- 0
    tmb_dat$S_p <- as(matrix(0, 1, 1), "sparseMatrix") 
    tmb_dat$S_p_n <- -1 
    map <- c(map, list(z_p = as.factor(NA), log_lambda_p = as.factor(NA)))
  } else {
    z_p <- rep(0, ncol(tmb_dat$Z_p)) 
    log_lambda_p <- rep(0, length(tmb_dat$S_p_n)) 
    random <- c(random, "z_p")
  }
  if (length(map) < 1) map <- NULL
  
  ## SETUP PARAMETERS FOR TMB  
  tmb_par <- list(beta_psi = beta_psi,
                  beta_p = beta_p, 
                  z_psi = z_psi, 
                  z_p = z_p, 
                  log_lambda_psi = log_lambda_psi, 
                  log_lambda_p = log_lambda_p)
  
  ## CREATE MODEL OBJECT 
  oo <- MakeADFun(data = tmb_dat, 
                  parameters = tmb_par, 
                  map = map, 
                  random = random,
                  DLL = "occu_tmb",
                  silent = !print)

  ## NORMALIZE  
  oo <- normalize(oo, flag = "flag")
  
  ## FIT MODEL
  fit <- nlminb(start = oo$par, objective = oo$fn, gradient = oo$gr)
  
  ## GET INFERENCE
  res <- sdreport(oo, getJointPrecision = TRUE)
  
  ## RESULTS
  val <- list(res = res, fit = fit, forms = forms, mats = mats)
  class(val) <- "occu"
  
  return(val)

}

#' Summary of occu object 
#'
#' @param obj output from fit_occu
#'
#' @return prints a summary for fixed effect parameters
#' @export
summary.occu <- function(obj) {
  return(print(obj$res))
}

#' Get predicted values for psi and p from parameters and matrices
#'
#' @param fix fixed parameter vector 
#' @param ran random effects parameter vector 
#' @param mats matrices from make_matrices
#'
#' @return list of two vectors: predicted values for psi and p 
get_predicted_values <- function(fix, ran, mats) {
  nms_fix <- names(fix)
  nms_ran <- names(ran)
  psi_pred <- mats$X_psi %*% fix[nms_fix == "beta_psi"] 
  if (!is.null(mats$Z_psi)) psi_pred <- psi_pred + mats$Z_psi %*% ran[nms_ran == "z_psi"]
  p_pred <- mats$X_p %*% fix[nms_fix == "beta_p"] 
  if (!is.null(mats$Z_p)) p_pred <- p_pred + mats$Z_p %*% ran[nms_ran == "z_p"]
  psi_pred <- plogis(psi_pred)
  p_pred <- plogis(p_pred)
  return(list(psi = psi_pred, p = p_pred))
}

#' Predict from fitted occupancy model
#'
#' @param obj fitted model object from fit_occu
#' @param visit_data see fit_occu
#' @param site_data see fit_occu
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also 
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @importFrom mgcv rmvn
predict.occu <- function(obj, visit_data, site_data, nboot = 0) {
  mats <- make_matrices(obj$forms, visit_data, site_data)
  fix <- obj$res$par.fixed
  ran <- obj$res$par.random
  pred <- get_predicted_values(fix, ran, mats)
  if (nboot > 0.5) {
    Q <- obj$res$jointPrecision
    Q <- Q[!grepl("log_lambda_", colnames(Q)), 
           !grepl("log_lambda_", colnames(Q)), drop = FALSE]
    V <- solve(Q)
    param <- c(fix, ran)
    param <- param[!grepl("log_lambda", names(param))]
    boots <- rmvn(nboot, param, V)
    colnames(boots) <- names(param)
    nfix <- length(fix[!grepl("log_lambda", names(fix))])
    boots_psi <- matrix(0, nr = nboot, nc = length(pred$psi))
    boots_p <- matrix(0, nr = nboot, nc = length(pred$p))
    for (b in 1:nboot) {
      boot_res <- get_predicted_values(boots[b, 1:nfix], boots[b, -(1:nfix)], mats)
      boots_psi[b,] <- boot_res$psi
      boots_p[b,] <- boot_res$p
    }
    pred$psiboot <- boots_psi
    pred$pboot <- boots_p
  }
  return(pred)
}

#' Compute estimated degrees of freedom for a spline 
#'
#' @param X design matrix for spline 
#' @param S smoothing matrix for spline 
#' @param lambda smoothiing parameters 
#' @param Sn number of parameters per spline 
#' 
#' @description This code was written by David L Miller as part of the CTMCdive R package
#' and is based on Wood et al. (2017) p 212. 
#'
#' @return estimated degrees of freedom 
edf.occu <- function(X, S, lambda, Sn){
  if (length(lambda) == 0) return(0)

  # duplicate lambda enough times
  lambda <- rep(lambda, Sn)
  
  # calculate lambda*S
  Sbig <- S * lambda
  
  # calculate the hat matrix
  XtX <- t(X) %*% X
  Fi <- solve(XtX + Sbig)
  F <- Fi %*% XtX
  
  # return the trace
  return(sum(diag(F)))
}

#' Compute (estimated) degrees of freedom for fitted occupancy model
#'
#' @param obj fitted model object 
#' @param by default is FALSE, if TRUE then list of degrees of freedom is returned, one fixed parameters, smooth
#' for psi, and smooth for p. When FALSE, these are summed and total estimated degrees of freedom returned. 
#'
#' @return degrees of freedom (see "by" above)
#' @export
dof.occu <- function(obj, by = FALSE) {
  fix <- obj$res$par.fixed
  df <- length(fix[!grepl("log_lambda", names(fix))])
  psi_lambda <- exp(fix[names(fix) == "log_lambda_psi"])
  psi_edf <- edf.occu(obj$mats$Z_psi, obj$mats$S_psi$S, psi_lambda, obj$mats$S_psi$Scols)
  p_lambda <- exp(fix[names(fix) == "log_lambda_p"])
  p_edf <- edf.occu(obj$mats$Z_p, obj$mats$S_p$S, p_lambda, obj$mats$S_p$Scols)
  if (by) return(list(fix = df, psi = psi_edf, p = p_edf))
  return(df + psi_edf + p_edf)
}

#' Log-likelihood for fitted occupancy model 
#'
#' @param object fitted occupancy model from fit_occu
#' @param ... 
#'
#' @return log-likelihood with estimated degrees of freedom as attribute "df"
#' @export
logLik.occu <- function(object, ...) {
  val <- -object$fit$objective 
  attributes(val)$df <- dof.occu(object)
  return(val)
}







