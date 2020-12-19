library(parallel)
library(future)
source("../utils/my_parLapply.R")
source("../utils/sims_params.R")
source("../utils/glmnet_sim_functions.R")

n_cores = future::availableCores()
cl = makeCluster(n_cores)
a=clusterCall(cl, function() library(glmnet))
a=clusterCall(cl, function() source("../utils/getdata.R"))
a=clusterCall(cl, function() source("../utils/glmnet_sim_functions.R"))
clusterSetRNGStream(cl, iseed=123)
param = param_p500k
alphas = c(1, 0.8, 0.6, 0.5, 0.4, 0.2, 0)
for (alpha in alphas) {
	for (i in 1:nrow(param)) {
		p = param$p[i]
		n = param$n[i]
		Nsim = param$Nsim[i]
		repl = glmnet_my_parLapply(cl, Nsim, n_cores, n, p, alpha)
		results = list()
		results$n = n
		results$p = p
		results$accuracies = simplify2array(repl)
		filename = paste0("sim_p", p, "_n", n, "_alpha", alpha, ".RDS")
		saveRDS(results, filename)
		print_summary(results)
	}
}

stopCluster(cl)