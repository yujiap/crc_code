library(here)
library(crc)
library(future)
library(parallel)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(glmnet)
library(pamr)
library(MASS)
library(sn)
source("sim_functions.R")

experiment_crc = function(params) {
  n = params$n
  p = params$p
  alpha.rs.level = params$alpha.rs.level
  L.type = params$L.type
  eps.type = params$eps.type
  min.eig.ratio = 0
  data = get_data(n, p, alpha.rs.level = alpha.rs.level, L.type = L.type, eps.type = eps.type)
  
  model = "simple"
  t.simple = system.time({crc.fit = crc(data$S, data$T, min.eig.ratio = min.eig.ratio)})['elapsed']
  acc.simple.S = get_theoretical_accuracy(data, model, crc.fit$beta1.S, crc.fit$beta0.S)
  acc.simple.L = get_theoretical_accuracy(data, model, crc.fit$beta1.L, crc.fit$beta0.L)
  acc.simple.C = get_theoretical_accuracy(data, model, crc.fit$beta1.C, crc.fit$beta0.C)
  nsel.simple = sum(crc.fit$S.fit$beta1.S.S != 0)
  sel.simple = which(crc.fit$S.fit$beta1.S.S != 0)
  
  model = "uncorrelated"
  t.uncor = system.time({crc.fit = crc(data$Z1, data$T, min.eig.ratio = min.eig.ratio)})['elapsed']
  acc.uncor.S = get_theoretical_accuracy(data, model, crc.fit$beta1.S, crc.fit$beta0.S)
  acc.uncor.L = get_theoretical_accuracy(data, model, crc.fit$beta1.L, crc.fit$beta0.L)
  acc.uncor.C = get_theoretical_accuracy(data, model, crc.fit$beta1.C, crc.fit$beta0.C)
  nsel.uncor = sum(crc.fit$S.fit$beta1.S.S != 0)
  sel.uncor = which(crc.fit$S.fit$beta1.S.S != 0)
  
  model = "correlated"
  t.corr = system.time({crc.fit = crc(data$Z2, data$T, min.eig.ratio = min.eig.ratio)})['elapsed']
  acc.corr.S = get_theoretical_accuracy(data, model, crc.fit$beta1.S, crc.fit$beta0.S)
  acc.corr.L = get_theoretical_accuracy(data, model, crc.fit$beta1.L, crc.fit$beta0.L)
  acc.corr.C = get_theoretical_accuracy(data, model, crc.fit$beta1.C, crc.fit$beta0.C)
  nsel.corr = sum(crc.fit$S.fit$beta1.S.S != 0)
  sel.corr = which(crc.fit$S.fit$beta1.S.S != 0)
  
  return(list(acc.simple.S=acc.simple.S, acc.simple.L=acc.simple.L, acc.simple.C=acc.simple.C,
              acc.uncor.S=acc.uncor.S, acc.uncor.L=acc.uncor.L, acc.uncor.C=acc.uncor.C,
              acc.corr.S=acc.corr.S, acc.corr.L=acc.corr.L, acc.corr.C=acc.corr.C,
              t.simple=t.simple, t.uncor=t.uncor, t.corr=t.corr,
              nsel.simple=nsel.simple, nsel.uncor=nsel.uncor, nsel.corr=nsel.corr,
              sel.simple=nsel.simple, sel.uncor=nsel.uncor, sel.corr=nsel.corr))
}

experiment_other = function(params) {
  n = params$n
  p = params$p
  alpha.rs.level = params$alpha.rs.level
  L.type = params$L.type
  eps.type = params$eps.type
  data = get_data(n, p, alpha.rs.level = alpha.rs.level, L.type = L.type, eps.type = eps.type)
  clf = match.fun(paste0("wrapper_", params$clf))
  
  model = "simple"
  t.simple = system.time({fit = clf(data$S, data$T)})['elapsed']
  acc.simple = get_theoretical_accuracy(data, model, fit$beta1, fit$beta0)
  nsel.simple = fit$n.selected
  sel.simple = fit$selected.ix
  if (params$clf == "glmnet") glmnet.lambda.simple = fit$lambda else glmnet.lambda.simple = NA
  if (params$clf == "pam") pam.best.thresh.simple = fit$best.thresh else pam.best.thresh.simple = NA
  
  model = "uncorrelated"
  t.uncor = system.time({fit = clf(data$Z1, data$T)})['elapsed']
  acc.uncor = get_theoretical_accuracy(data, model, fit$beta1, fit$beta0)
  nsel.uncor = fit$n.selected
  sel.uncor = fit$selected.ix
  if (params$clf == "glmnet") glmnet.lambda.uncor = fit$lambda else glmnet.lambda.uncor = NA
  if (params$clf == "pam") pam.best.thresh.uncor = fit$best.thresh else pam.best.thresh.uncor = NA
  
  model = "correlated"
  t.corr = system.time({fit = clf(data$Z2, data$T)})['elapsed']
  acc.corr = get_theoretical_accuracy(data, model, fit$beta1, fit$beta0)
  nsel.corr = fit$n.selected
  sel.corr = fit$selected.ix
  if (params$clf == "glmnet") glmnet.lambda.corr = fit$lambda else glmnet.lambda.corr = NA
  if (params$clf == "pam") pam.best.thresh.corr = fit$best.thresh else pam.best.thresh.corr = NA
  
  return(list(acc.simple=acc.simple, acc.uncor=acc.uncor, acc.corr=acc.corr, 
              t.simple=t.simple, t.uncor=t.uncor, t.corr=t.corr,
              nsel.simple=nsel.simple, nsel.uncor=nsel.uncor, nsel.corr=nsel.corr,
              sel.simple=sel.simple, sel.uncor=sel.uncor, sel.corr=sel.corr,
              glmnet.lambda.simple=glmnet.lambda.simple, pam.best.thresh.simple=pam.best.thresh.simple,
              glmnet.lambda.uncor=glmnet.lambda.uncor, pam.best.thresh.uncor=pam.best.thresh.uncor,
              glmnet.lambda.corr=glmnet.lambda.corr, pam.best.thresh.corr=pam.best.thresh.corr))
}

## Setup parallel exec
n_cores = future::availableCores()
cl = makeCluster(n_cores)
packages = c("crc", "magrittr", "here", "MASS", "sn", "glmnet", "pamr")
a = clusterExport(cl, list("experiment_crc", "experiment_other", "myDlda", "packages"))
a = clusterCall(cl, function() lapply(packages, library, character.only = TRUE))
a = clusterCall(cl, function() source(here("sim_functions.R")))
clusterSetRNGStream(cl, iseed = 123)

## Vary all parameters
sim_all = crossing(n = c(50, 100, 200, 500, 1000),
                   p = 100000,
                   alpha.rs.level = c(0.0001, 0.001, 0.01, 0.1, 1), 
                   L.type = c("normal", "mixture", "lognormal", "uniform"), 
                   eps.type = c("normal", "dexp", "unif", "t"))
sim_all = sim_all %>%
  mutate(nsim = case_when(n == 50 ~ 100,
                          n == 100 ~ 100,
                          n == 200 ~ 100,
                          n == 500 ~ 30,
                          n == 1000 ~ 10))

## Setup simulation
sim = sim_all %>% filter(eps.type == "normal", L.type == "normal")
sim = sim %>% mutate(acc.simple=NA, acc.uncor=NA, acc.corr=NA, 
                     t.simple=NA, t.uncor=NA, t.corr=NA,
                     nsel.simple=NA, nsel.uncor=NA, nsel.corr=NA,
                     sel.simple=NA, sel.uncor=NA, sel.corr=NA,
                     glmnet.lambda.simple=NA, pam.best.thresh.simple=NA,
                     glmnet.lambda.uncor=NA, pam.best.thresh.uncor=NA,
                     glmnet.lambda.corr=NA, pam.best.thresh.corr=NA)
sim = crossing(sim, clf=c("glmnet", "pam", "myDlda"))
sim = sim %>% arrange(clf)
results = as.list(rep(NA, nrow(sim)))
for (i in seq_len(nrow(sim))) {
  params = sim[i,]
  nsim = params$nsim
  print(i)
  print(params)
  
  res = chunked_parLapply(nsim, experiment_other, params, n_cores, cl)
  
  sim$acc.simple[[i]] = list(sapply(res, function(x) x$acc.simple))
  sim$acc.uncor[[i]]  = list(sapply(res, function(x) x$acc.uncor))
  sim$acc.corr[[i]]	 = list(sapply(res, function(x) x$acc.corr))
  
  sim$t.simple[[i]] = list(sapply(res, function(x) x$t.simple))
  sim$t.uncor[[i]]  = list(sapply(res, function(x) x$t.uncor))
  sim$t.corr[[i]]   = list(sapply(res, function(x) x$t.corr))
  
  sim$nsel.simple[[i]] = list(sapply(res, function(x) x$nsel.simple))
  sim$nsel.uncor[[i]]  = list(sapply(res, function(x) x$nsel.uncor))
  sim$nsel.corr[[i]]   = list(sapply(res, function(x) x$nsel.corr))
  
  sim$sel.simple[[i]] = list(lapply(res, function(x) x$sel.simple))
  sim$sel.uncor[[i]]  = list(lapply(res, function(x) x$sel.uncor))
  sim$sel.corr[[i]]   = list(lapply(res, function(x) x$sel.corr))
  
  sim$glmnet.lambda.simple[[i]] = list(sapply(res, function(x) x$glmnet.lambda.simple))
  sim$glmnet.lambda.uncor[[i]]  = list(sapply(res, function(x) x$glmnet.lambda.uncor))
  sim$glmnet.lambda.corr[[i]]   = list(sapply(res, function(x) x$glmnet.lambda.corr))
  
  sim$pam.best.thresh.simple[[i]] = list(sapply(res, function(x) x$pam.best.thresh.simple))
  sim$pam.best.thresh.uncor[[i]]  = list(sapply(res, function(x) x$pam.best.thresh.uncor))
  sim$pam.best.thresh.corr[[i]]   = list(sapply(res, function(x) x$pam.best.thresh.corr))
  
  sim = as_tibble(sim)
  saveRDS(sim, "sim_alpha_other.RDS")
}

stopCluster(cl)