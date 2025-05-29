library(SimDesign)
library(lme4)
library(lavaan)
library(tidyverse)
library(psych)
library(clubSandwich)
library(plyr)
library(arm)

# Conditions

DESIGNFACTOR <- createDesign(
  # average cluster size
  ave_clus_size = 30,
  # number of clusters
  num_clus = c(15, 30, 50),
  # icc
  icc = 0.1,
  # grand intercept
  gamma00 = 1,
  # gamma 10, within-cluster slope.
  gamma10 = c(0, 0.3),
  # gamma 01, slope for between-cluster predictor
  gamma01 = c(0, 0.3),
  # gamma 11, slope for cross-level interaction
  gamma11 = 0,
  # tau00sq = 1,
  balanced = TRUE,
  # variance pattern, 1 = homo, 2 = positive pairing, 3 = negative pairing
  lev1_vp = c(1, 2, 3),
  lev2_vp = c(1, 2, 3)
)
# add condition number
DESIGNFACTOR <- rowid_to_column(DESIGNFACTOR,"cond")

# Data Generations

Generate <- function(condition, fixed_objects = NULL) {
  # gamma 00 grand intercept
  gamma00 <- condition$gamma00
  # gamma 10, within-cluster slope.
  gamma10 <- condition$gamma10
  # gamma 01, slope for between-cluster predictor
  gamma01 <- condition$gamma01
  # gamma 11, slope for cross-level interaction
  gamma11 <- condition$gamma11
  # tau0sq <- condition$tau0sq
  ave_clus_size <- condition$ave_clus_size
  num_clus <- condition$num_clus
  icc <- condition$icc
  # total number of observation
  num_obs <- ave_clus_size * num_clus
  sigma = 1
  ## icc 
  if(icc == 0.1){
    tau0sq = (1/9) * sigma^2
    tau0 = sqrt(tau0sq)
  }else{
    tau0sq = (3/7) * sigma^2
    tau0 = sqrt(tau0sq)
  }
  tau1 = tau0
  tau1sq = tau1 ^ 2
  ## unbalanced design
  if (condition$balanced) {
    clus_size <- rep(ave_clus_size, 5)
  }else {
    clus_size <- seq(ave_clus_size * 2 / 6, 
                     ave_clus_size * 10 / 6, length.out = 5)
  }
  clus_id <- rep.int(seq_len(num_clus), rep_len(clus_size, num_clus))
  # Generated Level 2 predictor
  Z <- rnorm(num_clus)
  # expand to level 1
  Z <- Z[clus_id] 
  # Generated level 1 predictor
  X <- rnorm(num_obs)
  u0j <- rnorm(num_clus, sd = tau0)[clus_id]
  u1j <- rnorm(num_clus, sd = tau1)[clus_id]
  # simulate heteroscedasticity
  if(condition$lev2_vp == 1){
    lambdaz <- 1
  }else if(condition$lev2_vp == 2 && ave_clus_size == 6){
    lambdaz <- (abs(Z) + 1)* 0.23 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev2_vp == 2 && ave_clus_size == 30){
    lambdaz <- (abs(Z) + 1)* 0.10 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev2_vp == 3 && ave_clus_size == 6){
    lambdaz <- (1/(abs(Z) + 1))* 0.73 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev2_vp == 3 && ave_clus_size == 30){
    lambdaz <- (1/(abs(Z) + 1))* 0.30 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }
  eij <- rnorm(num_obs, sd = sigma) 
  if(condition$lev1_vp == 1){ 
    lambdax <- 1
  }else if(condition$lev1_vp == 2 && ave_clus_size == 6){
    lambdax <- (abs(X) + 1) * 0.23 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev1_vp == 2 && ave_clus_size == 30){
    lambdax <- (abs(X) + 1) * 0.10 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev1_vp == 3 && ave_clus_size == 6){
    lambdax <- (1/(abs(X) + 1)) * 0.73 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }else if(condition$lev1_vp == 3 && ave_clus_size == 30){
    lambdax <- (1/(abs(X) + 1)) * 0.30 * sqrt(rep_len(clus_size, num_clus)[clus_id])
  }
  
  # Generate observed score
  y <- gamma00 + gamma01*Z + gamma10 * X + gamma11*Z*X + lambdaz * u0j + lambdaz * u1j * X + lambdax * eij
  dat <- data.frame(y = y, z = Z,x = X, u0j = u0j, u1j = u1j, eij = eij, lev2_vi = lambdaz*u0j, lev2_vs = lambdaz*u1j*X, lev1_var = lambdax * eij, cid = clus_id)
  dat
}
# test
# dat1 <- Generate(DESIGNFACTOR[1,])
# dat2 <- Generate(DESIGNFACTOR[2,])
# dat100 <- Generate(DESIGNFACTOR[3,])

# Helper functions for fitting random slopes model with REML
run_rsREML <- function(data) {
  fit_rs <- lmer(y ~ z * x+ (x| cid), data = data)
  sumreml <- as.data.frame(coef(summary(fit_rs)))
  est <- sumreml["Estimate"]
  se <- sumreml["Std. Error"]
  tstat <- sumreml["t value"]
  pval <- sumreml["Pr(>|t|)"]
  # df <-sumKR["df"]
  c(est["(Intercept)",], se["(Intercept)",], est["z",],se["z",],tstat["z",], 
    pval["z",], est["x",],se["x",],tstat["x",], pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
# run_rsREML(dat1)

# Helper functions for fitting random intercepts model with KR
run_riKR <- function(data) {
  fit_ri <- lmer(y ~ z * x+ (1| cid), data = data)
  sumKR <- as.data.frame(coef(summary(fit_ri, ddf = "Kenward-Roger")))
  est <- sumKR["Estimate"]
  se <- sumKR["Std. Error"]
  tstat <- sumKR["t value"]
  pval <- sumKR["Pr(>|t|)"]
  # df <-sumKR["df"]
  c(est["(Intercept)",], se["(Intercept)",], est["z",],se["z",],tstat["z",], 
    pval["z",], est["x",],se["x",],tstat["x",], pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
run_riKR(dat1)

# Helper functions for fitting random slopes model with KR
run_rsKR <- function(data) {
  fit_rs <- lmer(y ~ z * x+ (x| cid), data = data)
  sumKR <- as.data.frame(coef(summary(fit_rs, ddf = "Kenward-Roger")))
  est <- sumKR["Estimate"]
  se <- sumKR["Std. Error"]
  tstat <- sumKR["t value"]
  pval <- sumKR["Pr(>|t|)"]
  # df <-sumKR["df"]
  c(est["(Intercept)",], se["(Intercept)",], est["z",],se["z",],tstat["z",], 
    pval["z",], est["x",],se["x",],tstat["x",], pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
# test
# run_rsKR(dat1)

# Helper functions for fitting CR2 to random intercept
run_ricr <- function(data){
  fit_ri <- lmer(y ~ z * x + (1| cid), data = data)
  sandt <- coef_test(fit_ri, vcov = "CR2", cluster = data$cid, test = "Satterthwaite")
  # convention <- coef_test(fit_r, vcov = "CR1", cluster = data$cid, test = "z")
  sandt <- as.data.frame(sandt)
  est <- sandt["beta"]
  se <- sandt["SE"]
  tstat <- sandt["tstat"]
  # df <- sandt["df"]
  pval <- sandt["p_Satt"]
  c(est["(Intercept)",], se["(Intercept)",],
    est["z",],se["z",],tstat["z",], pval["z",], 
    est["x",],se["x",],tstat["x",],pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
# test
# run_rsKR(dat1)

# Helper functions for fitting CR2 to random slope
run_rscr <- function(data){
  fit_rs <- lmer(y ~ z * x + (x| cid), data = data)
  sandt <- coef_test(fit_rs, vcov = "CR2", cluster = data$cid, test = "Satterthwaite")
  # convention <- coef_test(fit_r, vcov = "CR1", cluster = data$cid, test = "z")
  sandt <- as.data.frame(sandt)
  est <- sandt["beta"]
  se <- sandt["SE"]
  tstat <- sandt["tstat"]
  # df <- sandt["df"]
  pval <- sandt["p_Satt"]
  c(est["(Intercept)",], se["(Intercept)",],
    est["z",],se["z",],tstat["z",], pval["z",], 
    est["x",],se["x",],tstat["x",],pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
# test
# run_rscr(dat1)

# Helper functions for fitting CR2 to ols
run_olscr <- function(data){
  fit_o <- lm(y ~ z * x, data = data)
  sandt <- coef_test(fit_o, vcov = "CR2", cluster = data$cid, test = "Satterthwaite")
  # convention <- coef_test(fit_r, vcov = "CR1", cluster = data$cid, test = "z")
  sandt <- as.data.frame(sandt)
  est <- sandt["beta"]
  se <- sandt["SE"]
  tstat <- sandt["tstat"]
  # df <- sandt["df"]
  pval <- sandt["p_Satt"]
  c(est["(Intercept)",], se["(Intercept)",],
    est["z",],se["z",],tstat["z",], pval["z",], 
    est["x",],se["x",],tstat["x",],pval["x",],
    est["z:x",],se["z:x",],tstat["z:x",],pval["z:x",])
}
# test
# run_olscr(dat1)

# Analyze function
Analyse <- function(condition, dat, fixed_objects = NULL) {
  # Initialize output
  methods <- c("REML","KR RI", "KR RS", "Adjusted CRVE RI", "Adjusted CRVE RS", "Adjusted CRVE OLS")
  stats <- c("est", "SE","test_stat","p_val")
  pars <- c("gamma00", "gamma01", "gamma10", "gamma11")
  # condition <- c("Lev1_out", "Lev2_out","heteroscedasticity")
  ret <- vector("list",length(methods))
  names(ret) <- methods
  ret$REML <- run_rsREML(dat)
  ret$kr_ri <- run_riKR(dat)
  ret$kr_rs <- run_rsKR(dat)
  ret$adjusted_ri <- run_ricr(dat)
  ret$adjusted_rs <- run_rscr(dat)
  ret$adjusted_ols <- run_olscr(dat)
  # Reformat the output
  ret <- unlist(ret)
  names(ret) <- rbind(outer(as.vector(outer(stats[c(1,2)], pars[1], paste, sep = "_")), 
                            methods, paste, sep = "_"),
                      outer(as.vector(outer(stats, pars[2:4], paste, sep = "_")), 
                            methods, paste, sep = "_"))
  ret
}

# test
# results <-data.frame(rbind(Analyse(DESIGNFACTOR[1,], dat),Analyse(DESIGNFACTOR[2,], dat2)))

tr_rsebias <- function(est_se, est, trim = .1) {
  est_se <- as.matrix(est_se)
  est <- as.matrix(est)
  est_se <- apply(est_se, 2, mean, trim = trim, na.rm = TRUE)
  if (trim == 0) {
    emp_sd <- apply(est, 2L, sd, na.rm = TRUE)
  } else {
    emp_sd <- apply(est, 2L, mad, na.rm = TRUE)
  }
  est_se / emp_sd - 1
}

# Evaluation

Evaluate <- function(condition, results, fixed_objects = NULL) {
  methods <- c("REML", "KR RI", "KR RS", "Adjusted CRVE RI", "Adjusted CRVE RS", "Adjusted CRVE OLS")
  stats <- c("est", "SE","test_stat","p_val")
  pars <- c("gamma00", "gamma01", "gamma10", "gamma11")
  pop <- c(condition$gamma00, condition$gamma01, condition$gamma10, condition$gamma11)
  results_est <- results[,grep("^est\\_", colnames(results))]
  results_se <- results[,grep("^(SE)_", colnames(results))]
  results_p <- results[,grep("^(p_val)_", colnames(results))]
  # type1err <- results_p < .05
  # lessthan0.05 <- EDR(results_p, alpha = 0.05)
  detrate_sig <- colMeans(results_p < 0.05)
  ret <- c(bias = tr_rsebias(results_se, results_est, trim = 0), detection_rate_sig = detrate_sig)
  ret
}
# test
# Evaluate(DESIGNFACTOR[1:2,], results)

res <- runSimulation(design = DESIGNFACTOR,
                     replications = 1000,
                     generate = Generate,
                     analyse = Analyse,
                     summarise = Evaluate,
                     ncores = 10,
                     packages = c("lme4", "lavaan","clubSandwich", "lmerTest","arm","tidyverse","lmeresampler", "CR2", "pbkrtest"),
                     filename = "sim7_trial3",
                     seed = rep(23506200, nrow(DESIGNFACTOR)),
                     #allow_na = TRUE,
                     parallel = TRUE,
                     save = TRUE,
                     save_results = TRUE
)
saveRDS(res, "sim7_trial3.RDS")