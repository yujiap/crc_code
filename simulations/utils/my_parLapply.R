# A version of parLapply, seems to help reduce memory overhead

# See: https://r.789695.n4.nabble.com/mclapply-memory-leak-td4711759.html
chunked_parLapply = function(cl, Nsim, ncores, n, p) {
	X = 1:Nsim
	nchunks = Nsim/ncores
	splits = lapply(1:nchunks, function(i) rep(i, ncores))
	if (Nsim %% floor(nchunks) != 0) {
		splits[[(nchunks+1)]] = rep( (floor(nchunks)+1), Nsim %% floor(nchunks))
	}
	i.list = split(X, factor(unlist(splits)))
	result.list = list()
	for (i in seq_along(i.list)) {
		i.vec = i.list[[i]]
		result.list[i.vec] = parLapply(cl, X[i.vec], 
			function(n, p, i_sim) {
				sim = getdata(n,p)
	    			getaccuracies(sim)
			}, n=n, p=p)
	}
	return(result.list)
}

my_parLapply = function(cl, Nsim, ncores, n, p) {
	if (Nsim > ncores) {
		repl = chunked_parLapply(cl, Nsim, ncores, n, p)
	} else {
		repl = parLapply(cl, 1:Nsim, 
			function(n, p, i_sim) {
				sim = getdata(n,p)
	    			getaccuracies(sim)
			}, n=n, p=p)
	}
	return(repl)
}

glmnet_chunked_parLapply = function(cl, Nsim, ncores, n, p, alpha) {
	X = 1:Nsim
	nchunks = Nsim/ncores
	splits = lapply(1:nchunks, function(i) rep(i, ncores))
	if (Nsim %% floor(nchunks) != 0) {
		splits[[(nchunks+1)]] = rep( (floor(nchunks)+1), Nsim %% floor(nchunks))
	}
	i.list = split(X, factor(unlist(splits)))
	result.list = list()
	for (i in seq_along(i.list)) {
		i.vec = i.list[[i]]
		result.list[i.vec] = parLapply(cl, X[i.vec], 
			function(alpha, n, p, i_sim) {
				sim = getdata(n,p)
	    			getaccuracies(sim, alpha)
			}, n=n, p=p, alpha=alpha)
	}
	return(result.list)
}

glmnet_my_parLapply = function(cl, Nsim, ncores, n, p, alpha) {
	if (Nsim > ncores) {
		repl = glmnet_chunked_parLapply(cl, Nsim, ncores, n, p, alpha)
	} else {
		repl = parLapply(cl, 1:Nsim, 
			function(alpha, n, p, i_sim) {
				sim = getdata(n,p)
	    			getaccuracies(sim, alpha)
			}, n=n, p=p, alpha=alpha)
	}
	return(repl)
}