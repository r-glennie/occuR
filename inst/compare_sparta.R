# Compare occuR and sparta 

library(sparta)
library(occuR)
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

# formulae 
forms <- list(psi ~ s(occasion, bs = "cs", k = 5),  
              p ~ 1)

# name formulae (needed for make_matrices)
names(forms) <- c(as.character(forms[[1]][[2]]),
                  as.character(forms[[2]][[2]]))

# get model matrices
mats <- make_matrices(forms, visit_data, site_data)

# true parameters 
beta_psi <- c(0.4)
beta_p <- c(0.7)
z_psi <- rnorm(4, 0, 1)
z_p <- NULL #rnorm(13, 0, 0.5)

logit_psi <- (mats$X_psi %*% c(beta_psi, z_psi))[attr(mats$X_psi, "index")]
logit_p <- (mats$X_p %*% c(beta_p, z_p))[attr(mats$X_p, "index")] 
psi <- plogis(logit_psi)
p <- plogis(logit_p)

plot(site_data$occasion, psi)

## Simulate survey 
occu <- rbinom(nrow(site_data), 1, psi)
occu2 <- rep(occu, each = nvisits)
visit_data$obs <- rbinom(totvisits, 1, p * occu2)

## Fit basic model 
system.time(m0 <- fit_occu(list(psi ~ s(occasion, bs = "cs", k = 5), p ~ 1), visit_data, site_data)) 
m0

# Fit model in sparta  ----------------------------------------------------
date <- Sys.Date() + visit_data$occasion * 365 + visit_data$visit

dat <- data.frame(site = visit_data$site, 
                  occasion =  date,  
                  taxa = ifelse(visit_data$obs == 1, "a", "b"))

system.time(
spar <- occDetModel(taxa = dat$taxa,
                         site = dat$site,
                         survey = dat$occasion,
                         species_list = c('a'),
                         write_results = FALSE,
                         n_iterations = 40000,
                         burnin = 20000,
                         n_chains = 3,
                         thinning = 3,
                         seed = 123))


plot(spar$a)

sparpsi <- spar$a$BUGSoutput$mean$psi.fs
mlepsi <- predict(m0, visit_data, data.table(site = 1, occasion = 1:5))$psi

plot(sparpsi, lwd = 1.5, type = "b", pch = 19, col = "red")
lines(mlepsi, lwd = 1.5, type = "b", pch = 19, col = "blue")
