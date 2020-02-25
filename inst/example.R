## Occupancy example 

library(occuTMB)

## Setup survey 
nsites <- 100 # number of sites
nocc <- 5 # number of occasions
nvisits <- matrix(10, nr = nsites, nc = nocc) # number of visits per site x occasion 

# set up data
site_data <- CJ(site = 1:nsites, occasion = 1:nocc)
site_data <- site_data[order(occasion),]
totvisits <- sum(nvisits)
visit_data <- data.table(site = rep(1, totvisits), occasion = 1, visit = 1, y = 1)

# fill out visit_data
k <- 1
for (occ in 1:nocc) {
  for (site in 1:nsites) {
    v <- 1 
    for (visit in 1:nvisits[site, occ]) {
      visit_data$site[k] <- site
      visit_data$occasion[k] <- occ 
      visit_data$visit[k] <- v 
      v <- v + 1
      k <- k + 1
    }
  }
}

# add covariates
visit_data$temp <- rnorm(nrow(visit_data))
site_data$hab <- rep(factor(sample(1:3, size = nsites, prob = c(0.2, 0.4, 0.4), replace = TRUE)), nocc)
site_data$cover <- rnorm(nrow(site_data))

# formulae 
forms <- list(psi ~ s(occasion, bs = "cs", k = 3), 
              p ~ s(temp, bs = "cs"))

# name formulae (needed for make_matrices)
names(forms) <- c(as.character(forms[[1]][[2]]),
                  as.character(forms[[2]][[2]]))

# get model matrices
mats <- make_matrices(forms, visit_data, site_data)

# true parameters 
beta_psi <- c(0.4)
beta_p <- c(-0.4)
z_psi <- c(0.5, -0.2)
z_p <- rnorm(9)

logit_psi <- mats$X_psi %*% beta_psi 
if (!is.null(z_psi)) logit_psi <- logit_psi + mats$Z_psi %*% z_psi
logit_p <- mats$X_p %*% beta_p 
if (!is.null(z_p)) logit_p <- logit_p + + mats$Z_p %*% z_p
psi <- plogis(logit_psi)
p <- plogis(logit_p)

## Simulate survey 
occu <- rbinom(nsites * nocc, 1, psi)
visit_data$y <- rbinom(totvisits, 1, p * occu[visit_data$site + nsites * (visit_data$occasion - 1)])

## Fit model
mod <- fit_occu(forms, visit_data, site_data) 

## Look at estimates
summary(mod)

## Predict phi and p 
pred <- predict(mod, visit_data, site_data)

## Predict for certain covariates
new_site_data <- data.table(site = 1:5, occasion = 1:5)
pred2 <- predict(mod, visit_data, new_site_data, nboot = 100)

plot(1:5, pred2$psi, type = "b")
matlines(t(pred2$psiboot), col = "grey80")
ucl <- apply(pred2$psiboot, 2, quantile, probs = 0.975)
lcl <- apply(pred2$psiboot, 2, quantile, probs = 0.025)
plot(1:5, pred2$psi, type = "b", pch = 20, ylim = c(min(lcl), max(ucl)))
lines(1:5, lcl, lty = "dotted")
lines(1:5, ucl, lty = "dotted")

tempgrid <- seq(min(visit_data$temp), max(visit_data$temp), length = 100)
new_visit_data <- data.table(site = 1, occasion = 1, visit = 1, temp = tempgrid)

predtemp <- predict(mod, new_visit_data, site_data, nboot = 1000)
ucl <- apply(predtemp$pboot, 2, quantile, probs = 0.975)
lcl <- apply(predtemp$pboot, 2, quantile, probs = 0.025)
plot(tempgrid, predtemp$p, type = "l", lwd = 1.5, ylim = c(min(lcl), max(ucl)))
lines(tempgrid, lcl, lwd = 1.5, lty = "dotted")
lines(tempgrid, ucl, lwd = 1.5, lty = "dotted")

mats <- make_matrices(forms, new_visit_data, site_data)
truetemp <- plogis(mats$X_p %*% beta_p + mats$Z_p %*% z_p)
lines(tempgrid, truetemp, col = "red", lwd = 1.5)

## Test by simulation
nsims <- 100
mods <- vector(mode = "list", length = nsims)
for (sim in 1:nsims) {
  cat(sim, " / ", nsims, "\r")
  ## Simulate data 
  occu <- rbinom(nsites * nocc, 1, psi)
  visit_data$y <- rbinom(totvisits, 1, p * occu[visit_data$site + nsites * (visit_data$occasion - 1)])
  ## Fit model
  mods[[sim]] <- fit_occu(forms, visit_data, site_data, print = FALSE)
}

ests <- sapply(mods, FUN = function(x) {x$res$par.fixed})
rowMeans(ests)
zests <- sapply(mods, FUN = function(x) {x$res$par.random[names(x$res$par.random) == "z_psi"]})
rowMeans(zests)
