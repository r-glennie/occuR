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
#' @importFrom mgcv gam uniquecombs
#' @importFrom mgcv predict.gam
#' @importFrom methods as
make_matrices <- function(forms, visit_data, site_data) {
  ## occupancy model
  site_data$psi <- 1:nrow(site_data)
  gam_psi <- gam(forms[["psi"]], data = site_data, method = "REML")
  site_data <- site_data[, psi := NULL]
  # get design matrix
  Xfull <- predict(gam_psi, newdata = site_data, type = "lpmatrix")
  X_psi <- uniquecombs(Xfull)
  # get smoothing matrix
  S_psi <- make_smoothing_matrix(gam_psi)

  ## detection model
  visit_data$p <- 1:nrow(visit_data)
  gam_p <- gam(forms[["p"]], data = visit_data, method = "REML")
  visit_data <- visit_data[, p := NULL]
  # get design matrix
  Xfull <- predict(gam_p, newdata = visit_data, type = "lpmatrix")
  X_p <- uniquecombs(Xfull)
  # get smoothing matrix
  S_p <- make_smoothing_matrix(gam_p)

  ## results
  res <- list(X_psi = X_psi,
              S_psi = S_psi,
              nfix_psi = gam_psi$nsdf,
              X_p = X_p,
              S_p = S_p,
              nfix_p = gam_p$nsdf,
              gam_p = gam_p,
              gam_psi = gam_psi)

  return(res)
}

#' Fit occupancy model
#'
#' @param forms list of formulae for psi (occupancy) and p (detection)
#' @param visit_data data.table with row for each visit with a "obs" column for detection record
#' @param site_data data.table with row for each site x occasion
#' @param start named list of starting values for beta_psi and beta_p - fixed effects parameters
#' @param print if TRUE then possibly useful info is printed out
#' @return fitted occupied model object
#' @export
#' @importFrom data.table uniqueN
#' @importFrom TMB MakeADFun sdreport normalize
#' @useDynLib occu_tmb
fit_occu <- function(forms, visit_data, site_data, start = NULL, print = TRUE) {
  ## DATA
  # order data
  site_data <- site_data[order(site_data$site, site_data$occasion),]
  visit_data <- visit_data[order(visit_data$site, visit_data$occasion, visit_data$visit),]
  # check data
  check_inp(forms, visit_data, site_data, start, print)
  # number of sites
  nsites <- uniqueN(site_data$site)
  # number of occasions
  nocc <- site_data[, .(n = .N), .(site)]$n
  # total number of detections per site x occasion
  nvis <- visit_data[, .(totsite = sum(obs), nvisit = .N), .(site, occasion)]
  totsite <- nvis$totsite
  nvisit <- nvis$nvisit

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
                  y = visit_data$obs,
                  totsite = totsite,
                  nvisit = nvisit,
                  X_psi = mats$X_psi,
                  psi_ind = attr(mats$X_psi, "index") - 1,
                  S_psi = mats$S_psi$S,
                  S_psi_n = as.integer(mats$S_psi$Scols),
                  X_p = mats$X_p,
                  p_ind = attr(mats$X_p, "index") - 1,
                  S_p = mats$S_p$S,
                  S_p_n = as.integer(mats$S_p$Scols))

  ## PARAMETERS
  map <- list()
  random <- NULL
  # fixed effects
  beta_psi <- rep(0, mats$nfix_psi)
  beta_p <- rep(0, mats$nfix_p)
  if (!is.null(start)) {
    len <- min(length(start$beta_psi), beta_psi)
    beta_psi[1:len] <- start$beta_psi[1:len]
    len <- min(length(start$beta_p), beta_p)
    beta_p[1:len] <- start$beta_p[1:len]
  }

  # random effects
  # occupancy
  if (is.null(mats$S_psi)) {
    z_psi <- 0
    log_lambda_psi <- 0
    tmb_dat$S_psi <- as(matrix(0, 1, 1), "sparseMatrix")
    tmb_dat$S_psi_n <- -1
    map <- c(map, list(z_psi = as.factor(NA), log_lambda_psi = as.factor(NA)))
  } else {
    z_psi <- rep(0, ncol(mats$X_psi) - mats$nfix_psi)
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
    z_p <- rep(0, ncol(mats$X_p) - mats$nfix_p)
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
  class(val) <- "occuR"

  return(val)

}

#' Simple input checks for fit_occu
#'
#' @inheritParams fit_occu
#'
#' @return stops execution if input is invalid
check_inp <- function(forms, visit_data, site_data, start, print) {

  # check visit_data
  if (!("data.table" %in% class(visit_data))) stop("visit_data should be a data.table")
  if (!("obs" %in% colnames(visit_data))) stop("visit_data does not have a column named 'obs'")
  if (!("site" %in% colnames(visit_data))) stop("visit_data does not have a column named 'site'")
  if (!("occasion" %in% colnames(visit_data))) stop("visit_data does not have a column named 'occasion'")
  if (!("visit" %in% colnames(visit_data))) stop("visit_data does not have a column named 'visit'")
  if ("p" %in% colnames(visit_data)) stop("visit_data cannot have a column named 'p'")
  if (any(abs(visit_data$obs) > 1e-10 & abs(visit_data$obs - 1) > 1e-10)) stop("visit_data obs has entries that are not zero or one")

  # check site data
  num <- site_data[, .(max = max(occasion), n = uniqueN(occasion))]
  if (any(num$max != num$n)) stop("site_data has missing occasions or occasions are mis-numered")
  if ("psi" %in% colnames(site_data)) stop("site_data cannot have a column named 'psi'")

  # consistency between visit_data and site_data
  if (!all(visit_data$site %in% site_data$site)) stop("visit_data has sites not included in site_data")
  siteocc <- paste0(visit_data$site, visit_data$occasion)
  siteocc2 <- paste0(site_data$site, site_data$occasion)
  if (!all(siteocc %in% siteocc2)) stop("visit_data has sites with occasions not included in site_data")
  if (!all(siteocc2 %in% siteocc)) stop("site_data has sites with occasions where no visits are recorded in visit_data")

  # check start
  if (!is.null(start)) {
    if (!("beta_p" %in% names(start))) stop("start must have list entry named 'beta_p'")
    if (!("beta_psi" %in% names(start))) stop("start must have a list entry named 'beta_psi'")
    if (!is.numeric(start$beta_p) | !is.numeric(start$beta_psi)) stop("start vectors are invalid")
  }

  if (!is.logical(print)) stop("print must be TRUE or FALSE")

}

#' Summary of occuR object
#'
#' @param obj output from fit_occu
#'
#' @return prints a summary
#' @export
summary.occuR <- function(obj, conf.level = 0.95) {
  est <- obj$res$par.fixed
  nms <- names(est)
  sd <- sqrt(diag(obj$res$cov.fixed))
  alpha <- 1 - (1 - conf.level) / 2
  lcl <- est - qnorm(alpha) * sd
  ucl <- est + qnorm(alpha) * sd
  val <- data.frame(estimate = est, sd = sd, lcl = lcl, ucl = ucl)
  val <- val[1:(obj$mats$nfix_psi + obj$mats$nfix_p),]
  nms <- nms[1:(obj$mats$nfix_psi + obj$mats$nfix_p)]
  val <- signif(val, 4)
  nms[nms == "beta_psi"] <- paste("psi:", colnames(obj$mats$X_psi[,1:obj$mats$nfix_psi, drop=FALSE]))
  nms[nms == "beta_p"] <- paste("p:", colnames(obj$mats$X_p[,1:obj$mats$nfix_p, drop = FALSE]))
  rownames(val) <- nms
  print(val)
  dof <- round(as.numeric(dof.occuR(obj, each = TRUE)), 1)
  cat(paste0("Total edf: ", sum(dof), "\t smooth_psi: ", dof[2], "\t smooth_p: ", dof[3]))
  invisible(obj)
}

#' Print of occuR object
#'
#' @param obj output from fit_occu
#'
#' @return prints a summary
#' @export
print.occuR <- function(obj) {
  summary(obj)
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
  psi_pred <- (mats$X_psi %*% c(fix[nms_fix == "beta_psi"], ran[nms_ran == "z_psi"]))
  p_pred <- (mats$X_p %*% c(fix[nms_fix == "beta_p"], ran[nms_ran == "z_p"]))
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
#' @importFrom Matrix solve
predict.occuR <- function(obj, visit_data, site_data, nboot = 0) {
  mats <- obj$mats
  site_data$psi <- 1:nrow(site_data)
  mats$X_psi <- predict(mats$gam_psi, newdata = site_data, type = "lpmatrix")
  site_data <- site_data[, psi := NULL]
  visit_data$p <- 1:nrow(visit_data)
  mats$X_p <- predict(mats$gam_p, newdata = visit_data, type = "lpmatrix")
  visit_data <- visit_data[, p := NULL]
  fix <- obj$res$par.fixed
  ran <- obj$res$par.random
  pred <- get_predicted_values(fix, ran, mats)
  if (nboot > 0.5) {
    Q <- obj$res$jointPrecision
    if (!is.null(Q)) {
      Q <- Q[!grepl("log_lambda_", colnames(Q)),
             !grepl("log_lambda_", colnames(Q)), drop = FALSE]
      V <- solve(Q)
    } else {
      V <- obj$res$cov.fixed
    }
    Q <- Q[!grepl("log_lambda_", colnames(Q)),
           !grepl("log_lambda_", colnames(Q)), drop = FALSE]
    V <- solve(Q)
>>>>>>> 92a3f93 (predict and edf use Matrix::solve)
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
#' @importFrom Matrix solve diag
edf.occuR <- function(X, S, lambda, Sn){
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
#' @param each default is FALSE, if TRUE then list of degrees of freedom is returned, one fixed parameters, smooth
#' for psi, and smooth for p. When FALSE, these are summed and total estimated degrees of freedom returned.
#'
#' @return degrees of freedom (see "each" above)
#' @export
dof.occuR <- function(obj, each = FALSE) {
  fix <- obj$res$par.fixed
  df <- length(fix[!grepl("log_lambda", names(fix))])
  psi_lambda <- exp(fix[names(fix) == "log_lambda_psi"])
  psi_edf <- edf.occuR(obj$mats$X_psi[, -(1:obj$mats$nfix_psi)], obj$mats$S_psi$S, psi_lambda, obj$mats$S_psi$Scols)
  p_lambda <- exp(fix[names(fix) == "log_lambda_p"])
  p_edf <- edf.occuR(obj$mats$X_p[, -(1:obj$mats$nfix_p)], obj$mats$S_p$S, p_lambda, obj$mats$S_p$Scols)
  if (each) return(list(fix = df, psi = psi_edf, p = p_edf))
  return(df + psi_edf + p_edf)
}

#' Log-likelihood for fitted occupancy model
#'
#' @param object fitted occupancy model from fit_occu
#' @param ...
#'
#' @return log-likelihood with estimated degrees of freedom as attribute "df"
#' @export
logLik.occuR <- function(object, ...) {
  val <- -object$fit$objective
  attributes(val)$df <- dof.occuR(object)
  return(val)
}







