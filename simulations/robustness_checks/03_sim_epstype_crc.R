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
library(nimble)
source("sim_functions.R")

experiment_crc = function(params) {
  n = params$n
  ntst = params$ntst
  p = params$p
  alpha.rs.level = params$alpha.rs.level
  L.type = params$L.type
  eps.type = params$eps.type
  xp = unlist(params$xp, recursive = FALSE)
  data = get_data_trn_tst(n, ntst, p, alpha.rs.level = alpha.rs.level, L.type = L.type, eps.type = eps.type)
  
  model = "simple"
  t.simple = system.time({crc.fit = train_test_crc(data$S, data$T, data$Stst, data$Ttst, xp)})['elapsed']
  acc.simple.S = crc.fit$accS
  acc.simple.L = crc.fit$accL
  acc.simple.C = crc.fit$accC
  
  model = "uncorrelated"
  t.uncor = system.time({crc.fit = train_test_crc(data$Z1, data$T, data$Z1tst, data$Ttst, xp)})['elapsed']
  acc.uncor.S = crc.fit$accS
  acc.uncor.L = crc.fit$accL
  acc.uncor.C = crc.fit$accC
  
  model = "correlated"
  t.corr = system.time({crc.fit = train_test_crc(data$Z2, data$T, data$Z2tst, data$Ttst, xp)})['elapsed']
  acc.corr.S = crc.fit$accS
  acc.corr.L = crc.fit$accL
  acc.corr.C = crc.fit$accC
  
  return(list(acc.simple.S=acc.simple.S, acc.simple.L=acc.simple.L, acc.simple.C=acc.simple.C,
              acc.uncor.S=acc.uncor.S, acc.uncor.L=acc.uncor.L, acc.uncor.C=acc.uncor.C,
              acc.corr.S=acc.corr.S, acc.corr.L=acc.corr.L, acc.corr.C=acc.corr.C,
              t.simple=t.simple, t.uncor=t.uncor, t.corr=t.corr))
}

experiment_other = function(params) {
  n = params$n
  ntst = params$ntst
  p = params$p
  alpha.rs.level = params$alpha.rs.level
  L.type = params$L.type
  eps.type = params$eps.type
  xp = unlist(params$xp, recursive = FALSE)
  data = get_data_trn_tst(n, ntst, p, alpha.rs.level = alpha.rs.level, L.type = L.type, eps.type = eps.type)
  clf = match.fun(paste0("train_test_", params$clf))
  
  model = "simple"
  t.simple = system.time({fit = clf(data$S, data$T, data$Stst, data$Ttst, params=xp)})['elapsed']
  acc.simple = fit$acc
  
  model = "uncorrelated"
  t.uncor = system.time({fit = clf(data$Z1, data$T, data$Z1tst, data$Ttst, params=xp)})['elapsed']
  acc.uncor = fit$acc
  
  model = "correlated"
  t.corr = system.time({fit = clf(data$Z2, data$T, data$Z2tst, data$Ttst, params=xp)})['elapsed']
  acc.corr = fit$acc
  
  return(list(acc.simple=acc.simple, acc.uncor=acc.uncor, acc.corr=acc.corr, 
              t.simple=t.simple, t.uncor=t.uncor, t.corr=t.corr))
}

## Vary all parameters
sim_all = crossing(n = c(50, 100, 200, 500, 1000),
                   ntst = 500,
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

## Setup parallel exec
n_cores = future::availableCores()
cl = makeCluster(n_cores)
packages = c("crc", "magrittr", "here", "MASS", "sn", "glmnet", "pamr", "nimble")
a = clusterExport(cl, list("experiment_crc", "experiment_other", "myDlda", "packages"))
a = clusterCall(cl, function() lapply(packages, library, character.only = TRUE))
a = clusterCall(cl, function() source(here("sim_functions.R")))
clusterSetRNGStream(cl, iseed = 123)

## Setup simulation
sim = sim_all %>% filter(alpha.rs.level == 1, L.type == "normal")
sim = sim %>% mutate(acc.simple.S=NA, acc.simple.L=NA, acc.simple.C=NA,
                     acc.uncor.S=NA, acc.uncor.L=NA, acc.uncor.C=NA,
                     acc.corr.S=NA, acc.corr.L=NA, acc.corr.C=NA,
                     t.simple=NA, t.uncor=NA, t.corr=NA)
sim = sim %>% crossing(xp = list(list(min.eig.ratio = 0)))

results = as.list(rep(NA, nrow(sim)))
for (i in seq_len(nrow(sim))) {
  params = sim[i,]
  nsim = params$nsim
  print(i)
  print(params)

  res = chunked_parLapply(nsim, experiment_crc, params, n_cores, cl)

  sim$acc.simple.S[[i]] = list(sapply(res, function(x) x$acc.simple.S))
  sim$acc.simple.L[[i]] = list(sapply(res, function(x) x$acc.simple.L))
  sim$acc.simple.C[[i]] = list(sapply(res, function(x) x$acc.simple.C))
  sim$acc.uncor.S[[i]]  = list(sapply(res, function(x) x$acc.uncor.S))
  sim$acc.uncor.L[[i]]  = list(sapply(res, function(x) x$acc.uncor.L))
  sim$acc.uncor.C[[i]]  = list(sapply(res, function(x) x$acc.uncor.C))
  sim$acc.corr.S[[i]]   = list(sapply(res, function(x) x$acc.corr.S))
  sim$acc.corr.L[[i]]   = list(sapply(res, function(x) x$acc.corr.L))
  sim$acc.corr.C[[i]]   = list(sapply(res, function(x) x$acc.corr.C))

  sim$t.simple[[i]] = list(sapply(res, function(x) x$t.simple))
  sim$t.uncor[[i]]  = list(sapply(res, function(x) x$t.uncor))
  sim$t.corr[[i]]   = list(sapply(res, function(x) x$t.corr))

  sim = as_tibble(sim)
  saveRDS(sim, "sim_epstype.RDS")
}
stopCluster(cl)