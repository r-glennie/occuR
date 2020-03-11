## Occupancy example 

library(occuTMB)

## Setup survey 
nsites <- 100 # number of sites
nocc <- 5 # number of occasions
nvisits <- matrix(5, nr = nsites, nc = nocc) # number of visits per site x occasion 

# set up data
totvisits <- sum(nvisits)
visit_data <- CJ(site = 1:nsites, occasion = 1:nocc, visit = 1:max(nvisits))

visit_data <- visit_data[order(site, occasion, visit)]
site_data <- visit_data[, .(occasion = unique(occasion)), .(site)]

# add covariates
#visit_data$temp <- rnorm(nrow(visit_data))
#site_data$hab <- rep(factor(sample(1:3, size = nsites, prob = c(0.2, 0.4, 0.4), replace = TRUE)), nocc)
#site_data$cover <- rnorm(nrow(site_data))

site_locs <- data.frame(site = 1:nsites, x = rnorm(nsites), y = rnorm(nsites))
site_data$x <- site_locs$x[site_data$site]
site_data$y <- site_locs$y[site_data$site]

# formulae 
forms <- list(psi ~ 1,  
              p ~ s(occasion, bs = "ts", k = 5))

# name formulae (needed for make_matrices)
names(forms) <- c(as.character(forms[[1]][[2]]),
                  as.character(forms[[2]][[2]]))

# get model matrices
mats <- make_matrices(forms, visit_data, site_data)

# true parameters 
beta_psi <- c(0.4)
beta_p <- c(0.7)
z_psi <- NULL#rnorm(29, 0, 0.1)
z_p <- rnorm(4, 0, 1)

logit_psi <- (mats$X_psi %*% c(beta_psi, z_psi))[attr(mats$X_psi, "index")]
logit_p <- (mats$X_p %*% c(beta_p, z_p))[attr(mats$X_p, "index")] 
psi <- plogis(logit_psi)
p <- plogis(logit_p)

## Simulate survey 
occu <- rbinom(nrow(site_data), 1, psi)
occu2 <- rep(occu, as.vector(nvisits[nvisits != 0]))
visit_data$obs <- rbinom(totvisits, 1, p * occu2)

## Fit model
start <- list(beta_psi = beta_psi, beta_p = beta_p)
mod <- fit_occu(forms, visit_data, site_data, start = start) 

## Look at estimates
summary(mod)

## Predict phi and p
pred <- predict(mod, visit_data, site_data)

plot(p, pred$p)
abline(a = 0, b = 1)

plot(1:nocc, psi[seq(50, length(psi), by = 100)], type = "p", pch = 19)
lines(1:nocc, pred$psi[seq(50, length(psi), by = 100)], type = "b", pch = 19, col = "red")

## Predict for certain covariates
#new_site_data <- data.table(site = 1:10, occasion = 1:10)
#pred2 <- predict(mod, visit_data, new_site_data, nboot = 100)

#plot(1:10, pred2$psi, type = "b")

# matlines(t(pred2$psiboot), col = "grey80")
# ucl <- apply(pred2$psiboot, 2, quantile, probs = 0.975)
# lcl <- apply(pred2$psiboot, 2, quantile, probs = 0.025)
# plot(1:5, pred2$psi, type = "b", pch = 20, ylim = c(min(lcl), max(ucl)))
# lines(1:5, lcl, lty = "dotted")
# lines(1:5, ucl, lty = "dotted")
# 
# tempgrid <- seq(min(visit_data$temp), max(visit_data$temp), length = 100)
# new_visit_data <- data.table(site = 1, occasion = 1, visit = 1, temp = tempgrid)
# 
# predtemp <- predict(mod, new_visit_data, site_data, nboot = 1000)
# ucl <- apply(predtemp$pboot, 2, quantile, probs = 0.975)
# lcl <- apply(predtemp$pboot, 2, quantile, probs = 0.025)
# plot(tempgrid, predtemp$p, type = "l", lwd = 1.5, ylim = c(min(lcl), max(ucl)))
# lines(tempgrid, lcl, lwd = 1.5, lty = "dotted")
# lines(tempgrid, ucl, lwd = 1.5, lty = "dotted")
# 
# mats <- make_matrices(forms, new_visit_data, site_data)
# truetemp <- plogis(mats$X_p %*% beta_p + mats$Z_p %*% z_p)
# lines(tempgrid, truetemp, col = "red", lwd = 1.5)
# 
# ## Test by simulation
nsims <- 100
mods <- vector(mode = "list", length = nsims)
for (sim in 1:nsims) {
  cat(sim, " / ", nsims, "\r")
  ## Simulate data
  occu <- rbinom(nrow(site_data), 1, psi)
  occu2 <- rep(occu, as.vector(nvisits[nvisits != 0]))
  visit_data$obs <- rbinom(totvisits, 1, p * occu2)
  ## Fit model
  mods[[sim]] <- fit_occu(forms, visit_data, site_data, print = FALSE)
}

ests <- sapply(mods, FUN = function(x) {x$res$par.fixed})
rowMeans(ests)
