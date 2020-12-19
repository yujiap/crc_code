library(parallel)
library(future)
source("../utils/my_parLapply.R")
source("../utils/sims_params.R")
source("../utils/pam_sim_functions.R")

n_cores = future::availableCores()
cl = makeCluster(n_cores)
a=clusterCall(cl, function() library(pamr))
a=clusterCall(cl, function() source("../utils/getdata.R"))
a=clusterCall(cl, function() source("../utils/pam_sim_functions.R"))
clusterSetRNGStream(cl, iseed=123)
param = param_p100k
for (i in 1:nrow(param)) {
	p = param$p[i]
	n = param$n[i]
	Nsim = param$Nsim[i]
	repl = my_parLapply(cl, Nsim, n_cores, n, p)
	results = list()
	results$n = n
	results$p = p
	results$accuracies = simplify2array(repl)
	filename = paste0("sim_p", p, "_n", n, ".RDS")
	saveRDS(results, filename)
	print_summary(results)
}

stopCluster(cl)