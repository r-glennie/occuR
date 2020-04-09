library(occuR)
library(unmarked)
library(ggplot2)

## load raw data (see Kery et al. 2013) 
data(crossbill)

## format data
# site data 
dat1 <- crossbill[apply(crossbill[,5:7], 1, FUN = function(x) {any(!is.na(x))}),]
site_data <- data.table(site = dat1$id, occasion = 1, elevation = dat1$ele, forest = dat1$forest)
for (i in 1:8) {
  dat1<- crossbill[apply(crossbill[,(5 + 3 * i):(7 + 3 *i)], 1, FUN = function(x) {any(!is.na(x))}),]
  site_data <- rbind(site_data, data.table(site = dat1$id, occasion = 1 + i, elevation = dat1$ele, forest = dat1$forest))
}
# scale all covariates
site_data$elevation <- scale(site_data$elevation)
site_data$forest <- scale(site_data$forest)
# visit data 
date <- as.matrix(crossbill[,32:58])
visdat <- as.matrix(crossbill[,5:31])
# set observations to NA if date covariate is missing 
visdat[is.na(date) != is.na(visdat)] <- NA
visit_data <- CJ(site = 1:uniqueN(site_data$site), occasion = 1:9, visit = 1:3)
visit_data$obs <- as.numeric(t(visdat))
# scale date covariate 
sddate <- sd(c(date), na.rm = TRUE)
meandate <- mean(date, na.rm = TRUE)
date <- (date - meandate) / sddate
visit_data$date <- as.numeric(t(date))
visit_data <- visit_data[!is.na(obs)]
# add forest/elevation covariates 
mch <- match(visit_data$site, site_data$site)
visit_data$forest <- site_data$forest[mch]
visit_data$elevation <- site_data$elevation[mch]


## fit models 

# fit basic model 
m0 <- fit_occu(list(psi ~ 1, p ~ 1), visit_data, site_data)

# DETECTION MODEL 
# forest effect 
mfor <- fit_occu(list(psi ~ 1, p ~ forest), visit_data, site_data)
mfors <- fit_occu(list(psi ~ 1, p ~ s(forest, bs = "cs")), visit_data, site_data)
mforquad <- fit_occu(list(psi ~ 1, p ~ forest + I(forest^2)), visit_data, site_data)

AIC(m0, mfor, mfors, mforquad)

forgr <- seq(-1.2, 2.3, 0.1)
forpred <- predict(mfors, data.table(site = 1, occasion = 1, visit = 1, date = 0, forest = forgr), site_data, nboot = 100)
ci <- apply(forpred$pboot, 2, quantile, prob = c(0.025, 0.975))
plot(forgr, forpred$p, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Forest cover", ylab = "Detection probability")
lines(forgr, ci[1,], lty = "dashed")
lines(forgr, ci[2,], lty = "dashed")

# elevation effect 
melev <- fit_occu(list(psi ~ 1, p ~ s(forest, bs = "cs") + s(elevation, bs = "cs")), visit_data, site_data)
melevquad <- fit_occu(list(psi ~ 1, p ~ s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data, 
                      start = list(beta_psi = 0.414, beta_p = -1.314))

AIC(mfors, melev, melevquad)

elevgr <- seq(-1.45, 2.43, 0.1)
elevpred <- predict(melevquad, data.table(site = 1, occasion = 1, visit = 1, date = 0, forest = 0, elevation = elevgr), site_data, nboot = 100)
ci <- apply(elevpred$pboot, 2, quantile, prob = c(0.025, 0.975))
plot(elevgr, elevpred$p, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Elevation", ylab = "Detection probability")
lines(elevgr, ci[1,], lty = "dashed")
lines(elevgr, ci[2,], lty = "dashed")

# year effect 
myr <- fit_occu(list(psi ~ 1, p ~ s(occasion, bs = "cs", k = 9) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)

myrquad <- fit_occu(list(psi ~ 1, p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)

AIC(myr, myrquad, melevquad)

yrgr <- seq(1, 9)
yrpred <- predict(myrquad, data.table(site = 1, occasion = yrgr, visit = 1, date = 0, forest = 0, elevation = 0), site_data, nboot = 100)
ci <- apply(yrpred$pboot, 2, quantile, prob = c(0.025, 0.975))
plot(yrgr + 1998, yrpred$p, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Year", ylab = "Detection probability")
lines(yrgr + 1998, ci[1,], lty = "dashed")
lines(yrgr + 1998, ci[2,], lty = "dashed")

# date effect 
mdate <- fit_occu(list(psi ~ 1, p ~ s(date, bs = "cs") + occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)

AIC(mdate, myrquad)

# OCCUPANCY MODEL

# forest effect 
mf <- fit_occu(list(psi ~ s(forest, bs = "cs"), p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)
mf2 <- fit_occu(list(psi ~ forest + I(forest^2), p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)

AIC(myrquad, mf, mf2)

forpred <- predict(mf, visit_data, data.table(site = 1, occasion = 1, forest = forgr), nboot = 100)
ci <- apply(forpred$psiboot, 2, quantile, prob = c(0.025, 0.975))
plot(forgr, forpred$psi, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Forest cover", ylab = "Occupancy probability")
lines(forgr, ci[1,], lty = "dashed")
lines(forgr, ci[2,], lty = "dashed")

# elevation effect 
me <- fit_occu(list(psi ~ s(forest, bs = "cs") + s(elevation, bs = "cs"), p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)
me2 <-  fit_occu(list(psi ~ s(forest, bs = "cs") + elevation + I(elevation^2), p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)
AIC(mf, me, me2)

elevpred <- predict(me2, visit_data, data.table(site = 1, occasion = 1,forest = 0, elevation = elevgr), nboot = 100)
ci <- apply(elevpred$psiboot, 2, quantile, prob = c(0.025, 0.975))
plot(elevgr, elevpred$psi, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Elevation", ylab = "Occupancy probability")
lines(elevgr, ci[1,], lty = "dashed")
lines(elevgr, ci[2,], lty = "dashed")

# year effect 
mt <-  fit_occu(list(psi ~ s(occasion, bs = "cs", k = 9) + s(forest, bs = "cs") + elevation + I(elevation^2), p ~ occasion + I(occasion^2) + s(forest, bs = "cs") + elevation + I(elevation^2)), visit_data, site_data)

AIC(mt, me2)

yrpred <- predict(mt, visit_data, data.table(site = 1, occasion = yrgr, date = 0, forest = 0, elevation = 0), nboot = 100)
ci <- apply(yrpred$psiboot, 2, quantile, prob = c(0.025, 0.975))
plot(yrgr + 1998, yrpred$psi, type = "l", lwd = 1.5, ylim = c(min(ci[1,]), max(ci[2,])), xlab = "Year", ylab = "Occupancy probability")
lines(yrgr + 1998, ci[1,], lty = "dashed")
lines(yrgr + 1998, ci[2,], lty = "dashed")
