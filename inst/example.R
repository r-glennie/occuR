## Occupancy example 

library(occuR)
library(ggplot2)
set.seed(53930)

## Setup survey 
nsites <- 100 # number of sites
nocc <- 5 # number of occasions
nvisits <- 10 # number of visits per site x occasion 

# set up data
totvisits <- nvisits * nsites * nocc
visit_data <- CJ(site = 1:nsites, occasion = 1:nocc, visit = 1:nvisits)
visit_data <- visit_data[order(site, occasion, visit)]
site_data <- visit_data[, .(occasion = unique(occasion)), .(site)]

# add covariates
visit_data$temp <- rnorm(nrow(visit_data), 12, 5)
site_data$hab <- rep(factor(sample(c("arable", "woodland", "urban"), size = nsites, prob = c(0.2, 0.4, 0.4), replace = TRUE)), nocc)
site_locs <- data.frame(site = 1:nsites, x = rnorm(nsites), y = rnorm(nsites))
site_data$x <- site_locs$x[site_data$site]
site_data$y <- site_locs$y[site_data$site]
site_hab <- site_data[occasion == 1]
visit_data$hab <- site_hab$hab[visit_data$site]

# formulae 
forms <- list(psi ~ t2(x, y, occasion, d = c(2, 1), bs = c("ts", "cs")) + hab,  
              p ~ s(occasion, bs = "cs", k = 5) + s(temp, bs = "cs"))

# name formulae (needed for make_matrices)
names(forms) <- c(as.character(forms[[1]][[2]]),
                  as.character(forms[[2]][[2]]))

# get model matrices
mats <- make_matrices(forms, visit_data, site_data)

# true parameters 
beta_psi <- c(0.4, 0.5, -1)
beta_p <- c(0.7)
z_psi <- rnorm(125, 0, 0.1)
z_p <- rnorm(13, 0, 0.5)

logit_psi <- (mats$X_psi %*% c(beta_psi, z_psi))[attr(mats$X_psi, "index")]
logit_p <- (mats$X_p %*% c(beta_p, z_p))[attr(mats$X_p, "index")] 
psi <- plogis(logit_psi)
p <- plogis(logit_p)

## Simulate survey 
occu <- rbinom(nrow(site_data), 1, psi)
occu2 <- rep(occu, each = nvisits)
visit_data$obs <- rbinom(totvisits, 1, p * occu2)

## Fit basic model 
m0 <- fit_occu(list(psi ~ 1, p ~ 1), visit_data, site_data) 
m0

## Habitat effect 
m_psihab <- fit_occu(list(psi ~ hab, p ~ 1), visit_data, site_data)
m_phab <- fit_occu(list(psi ~ 1, p ~ hab), visit_data, site_data)
m_psiphab <- fit_occu(list(psi ~ hab, p ~ hab), visit_data, site_data)
AIC(m0, m_psihab, m_phab, m_psiphab)

## Temperature effect 
m_temp <- fit_occu(list(psi ~ hab, p ~ temp), visit_data, site_data)
m_temps <- fit_occu(list(psi ~ hab, p ~ s(temp, bs = "cs")), visit_data, site_data)
AIC(m_psihab, m_temp, m_temps)

## Temporal effect 
m_pt <- fit_occu(list(psi ~ hab, p ~ s(temp, bs = "cs") + s(occasion, bs = "cs", k = 5)), visit_data, site_data)
m_psit <- fit_occu(list(psi ~ hab + s(occasion, bs = "cs", k = 5), p ~ s(temp, bs = "cs") + s(occasion, bs = "cs", k = 5)), visit_data, site_data)

AIC(m_temps, m_pt, m_psit)

## Spatio-Temporal effect 
m_psi_xyt <- fit_occu(list(psi ~ t2(x, y, occasion, bs = c("ts", "cs"), d = c(2, 1)) + hab, 
                           p ~ s(temp, bs = "cs") + s(occasion, bs = "cs", k = 5)), visit_data, site_data)

AIC(m_pt, m_psit, m_psi_xyt)

## Inference
m <- m_psi_xyt
m

# temperature effect 
tempgr <- seq(-5, 25, 0.1)
pred_temp <- predict(m, data.table(occasion = 1, temp = tempgr), site_data, nboot = 1000)
ci <- apply(pred_temp$pboot, 2, quantile, prob = c(0.025, 0.975))
plot(tempgr, pred_temp$p, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Temperature", ylab = "Detection Probability")
lines(tempgr, ci[1,], lty = "dotted")
lines(tempgr, ci[2,], lty = "dotted")

# occasion effect 
pred_occ <- predict(m, data.table(occasion = 1:nocc, temp = 18), site_data, nboot = 1000)
ci <- apply(pred_occ$pboot, 2, quantile, prob = c(0.025, 0.975))
plot(1:nocc, pred_occ$p, type = "b", pch = 19, lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Temperature", ylab = "Detection Probability")
lines(1:nocc, ci[1,], lty = "dotted")
lines(1:nocc, ci[2,], lty = "dotted")

# spatial effect 
xgr <- seq(-2, 2.5, 0.1)
ygr <- seq(-2, 2.5, 0.1)
gr <- expand.grid(xgr, ygr)
pred_xy <- predict(m, visit_data, data.table(occasion = 1, x = gr[,1], y = gr[,2], hab = "arable"), nboot = 1000)

ggplot() + 
  geom_tile(aes(x = gr[,1], y = gr[,2], fill = pred_xy$psi)) + 
  theme_bw() + 
  scale_x_continuous("x") + 
  scale_y_continuous("y") + 
  scale_fill_viridis_c("Occupancy")

# spatio-temporal effect 
xgr <- rep(gr[,1], nocc)
ygr <- rep(gr[,2], nocc)
tgr <- rep(1:nocc, each = nrow(gr))
pred_xyt <- predict(m, visit_data, data.table(occasion = tgr, x = xgr, y = ygr, hab = "arable"), nboot = 1000)

ggplot(data.frame(x = xgr, y = ygr, t = tgr, psi = pred_xyt$psi)) + 
  geom_tile(aes(x = x, y = y, group = t, fill = psi)) + 
  theme_bw() + 
  facet_wrap(~t) + 
  scale_x_continuous("x") + 
  scale_y_continuous("y") + 
  scale_fill_viridis_c("Occupancy")
