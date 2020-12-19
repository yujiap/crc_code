# Simulation functions modified to track which features selected

library(pamr)
library(glmnet)
library(crc)
library(magrittr)
library(parallel)

# --------------------------------------------------------------
# Wrappers for classifiers
# --------------------------------------------------------------

subgroupAcc = function(pred, Ttst, fctr_tst) {
	ret = rep(NA, length(levels(fctr_tst)))
	names(ret) = levels(fctr_tst)
	for (lvl in levels(fctr_tst)) {
		lvl.ix = which(fctr_tst==lvl)
		if (length(lvl.ix)!=0) {
			lvl.acc = mean((pred==Ttst)[lvl.ix])
			ret[lvl] = lvl.acc	
		}
	}
	return(ret)
}

crc_predict = function(fit, Ztst, clf=c("S", "L", "C")) {
  beta0 = fit[[paste0("beta0.", clf)]]
  beta1 = fit[[paste0("beta1.", clf)]]
  scores = Ztst %*% matrix(beta1) + beta0 
  pred = ifelse(scores < 0, fit$class_names[1], fit$class_names[2])
  return(pred)
}

getaccuracies.crc = function(Ztrn, T_trn, Ztst, Ttst, fctr_tst) {
  lvls = levels(fctr_tst)
  n.lvls = length(lvls)
  results = rep(NA, (4 + 3*n.lvls))
  names(results) = c("S", "L", "C", "n.selected", 
  					paste0("S.", lvls),
  					paste0("L.", lvls),
  					paste0("C.", lvls))
  fit = crc(Ztrn, T_trn, min.eig.ratio = 0) ###
  S_pred = crc_predict(fit, Ztst, "S")
  L_pred = crc_predict(fit, Ztst, "L")
  C_pred = crc_predict(fit, Ztst, "C")  
  S_acc = mean(S_pred==Ttst)
  L_acc = mean(L_pred==Ttst)
  C_acc = mean(C_pred==Ttst)
  results[1:3] = c(S_acc, L_acc, C_acc)
  results[4] = length(fit$S.selected)

  # Accuracies w/in levels of fctr
  results[paste0("S.", lvls)] = subgroupAcc(S_pred, Ttst, fctr_tst)
  results[paste0("L.", lvls)] = subgroupAcc(L_pred, Ttst, fctr_tst)
  results[paste0("C.", lvls)] = subgroupAcc(C_pred, Ttst, fctr_tst)

  ret = list(results = results, selected.ix = fit$S.selected)
  return(ret)
}


getaccuracies.glmnet = function(Ztrn, Ttrn, Ztst, Ttst, fctr_tst) {
  lvls = levels(fctr_tst)
  n.lvls = length(lvls)
  results = rep(NA, (3 + n.lvls))
  names(results) = c("glmnet", "lambda.1se", "n.selected", lvls)
  cv.fit = cv.glmnet(Ztrn, Ttrn, family="binomial", alpha=1, standardize=TRUE)
  pred = predict(cv.fit, newx=Ztst, type="class")
  results[1] = mean(pred == Ttst)
  results[2] = cv.fit$lambda.1se
  bestlam = which(cv.fit$lambda==cv.fit$lambda.1se)
  selected.ix = which(cv.fit$glmnet.fit$beta[,bestlam]!=0)
  results[3] = length(selected.ix)
  results[lvls] = subgroupAcc(pred, Ttst, fctr_tst)
  ret = list(results = results, selected.ix = selected.ix)
  return(ret)
}


getaccuracies.pam = function(Ztrn, Ttrn, Ztst, Ttst, fctr_tst) {
  lvls = levels(fctr_tst)
  n.lvls = length(lvls)
  results = rep(NA, (3 + n.lvls))
  names(results) = c("pam", "n.selected", "best.thresh", lvls)
  cv.fit = pamr.train(list("x"=t(Ztrn), "y"=Ttrn))
  t0 = which(cv.fit$errors == min(cv.fit$errors)) %>% max
  best.thresh = cv.fit$threshold[t0]
  pred = pamr.predict(cv.fit, newx=t(Ztst), threshold=best.thresh)
  results[1] = mean(pred == Ttst)
  results[2] = cv.fit$nonzero[t0]
  results[3] = best.thresh
  selected.ix = pamr.listgenes(fit=cv.fit, list("x"=t(Ztrn), "y"=Ttrn, "geneid"=1:ncol(Ztrn)), threshold=best.thresh)[,"id"]
  results[lvls] = subgroupAcc(pred, Ttst, fctr_tst)
  ret = list(results = results, selected.ix = selected.ix)
  return(ret)
}


getaccuracies.myDlda = function(Ztrn, Ttrn, Ztst, Ttst, fctr_tst) {
  lvls = levels(fctr_tst)
  n.lvls = length(lvls)
  results = rep(NA, (2 + n.lvls))
  names(results) = c("myDlda", "n.selected", lvls)
  fit = myDlda(Ztrn, Ttrn)
  pred = myDlda_predict(fit, Ztst)$pred
  results[1] = mean(pred == Ttst)
  results[2] = length(fit$selected.idx)
  results[lvls] = subgroupAcc(pred, Ttst, fctr_tst)
  ret = list(results = results, selected.ix = fit$selected.idx)
  return(ret)
}


myDlda = function(S, class, prior = c(0.5, 0.5), S.grid = crc:::get_S.grid("default", S)) {
  I = crc:::get_I(class)
  n = length(class)
  n1 = sum(I[,1])
  n2 = sum(I[,2])
  
  Cs = t(I) %*% S
  Cs.sq = t(I) %*% S^2
  Mu1 = Cs[1,] / n1
  Mu2 = Cs[2,] / n2
  ss1 = Cs.sq[1,] - Cs[1,]^2/n1
  ss2 = Cs.sq[2,] - Cs[2,]^2/n2
  
  nN = length(S.grid)
  scores.LOO_N = matrix(0, n, nN)
  for (i in 1:n) {
    sigmas2 = crc:::downdatePooledVar(S, I, Cs, Cs.sq, ss1, ss2, i)
    mu1 = if (I[i,1]) (Cs[1,] - S[i,])/(n1-1) else Mu1
    mu2 = if (I[i,2]) (Cs[2,] - S[i,])/(n2-1) else Mu2
    n1.n2 = colSums(I[-i,])
    tdrop = ((mu2 - mu1) / sqrt(sigmas2)) * sqrt((n1.n2[1]*n1.n2[2])/(n-1))
    toporder = order(abs(tdrop), decreasing = T)
    for (j in 1:nN) {
      N = S.grid[j]
      topN = toporder[1:N]
      beta1 = (mu2[topN] - mu1[topN]) / sigmas2[topN]
      beta0 = log(prior[2]/prior[1]) - 0.5 * sum((mu2[topN] + mu1[topN])*beta1)
      scores.LOO_N[i,j] = S[i, topN, drop = F] %*% beta1 + beta0
    }
  }
  
  thr_acc = apply(scores.LOO_N, 2,
                  function(x) {
                    sm1 = mean(x[I[,1] == 1])
                    sm2 = mean(x[I[,2] == 1])
                    spv = (n1*var(x[I[,1] == 1]) + n2*var(x[I[,2] == 1]))/(n-2)
                    p.mistake = prior[1]*pnorm(sm1/sqrt(spv)) + prior[2]*pnorm(-sm2/sqrt(spv))
                    acc = 1 - p.mistake
                  })
  
  sigmas2 = (ss1 + ss2) / (n-2)
  t = (Mu2 - Mu1) / sqrt(sigmas2 * n/(n1*n2))
  
  idx = which.max(thr_acc)
  topN = order(abs(t), decreasing = T)[1:S.grid[idx]]
  beta1 = (Mu2[topN] - Mu1[topN]) / sigmas2[topN]
  beta0 = log(prior[2]/prior[1]) - 0.5*sum((Mu2[topN] + Mu1[topN])*beta1)

  return(list(
              "beta1" = beta1,
              "beta0" = beta0,
              "scores.LOO"=scores.LOO_N[,idx],
              "selected.idx"=topN,
              "v"=thr_acc,
              "class_names" = colnames(I)
  ))
}

myDlda_predict = function(fit, Ytest) {
  scores = Ytest[, fit$selected.idx, drop = F] %*% fit$beta1 + fit$beta0
  pred = ifelse(scores > 0, fit$class_names[2], fit$class_names[1])
  return(list('scores'=scores, 'pred'=pred))
}


# --------------------------------------------------------------
# Simulation fn, for clustered data (e.g. multiple samples per id)
# --------------------------------------------------------------

# Train/test split for clustered observations which keeps clusters intact (i.e. no splitting clusters
# across train/test sets)
sample_clusters = function(id_clusters, id_class, train_size) {
	classes = factor(id_class)
	n1=sum(classes==levels(classes)[1])
	n2=sum(classes==levels(classes)[2])
	n0=min(n1,n2)
	
	if (train_size < 1) n0trn = floor(n0*train_size)
	if (train_size >= 1) n0trn = train_size
    if (n0trn >= n0) {print("train_size too large"); return()}
    n0tst = n0 - n0trn

	sample1 = sample(which(classes==levels(classes)[1]), n0)
	sample2 = sample(which(classes==levels(classes)[2]), n0)
	clusters_trn = c(sample1[1:n0trn], sample2[1:n0trn])
	clusters_tst = c(sample1[(n0trn+1):n0], sample2[(n0trn+1):n0])
	samples_trn = unlist(id_clusters[clusters_trn])
	samples_tst = unlist(id_clusters[clusters_tst])

	return(list("clusters_trn"=clusters_trn,
				"clusters_tst"=clusters_tst,
				"samples_trn"=samples_trn, 
				"samples_tst"=samples_tst,
				"n0_trn"=n0trn,
				"n0tst"=n0tst))
}

run_sim_clustered = function(clf, predictors, classes, id, fctr, Nsim, downsample, train_size, filename) {
	GETACCURACIES = match.fun(paste0("getaccuracies.", clf))
	wrappers = c("crc_predict",
			"myDlda",
			"myDlda_predict", 
			"getaccuracies.crc", 
			"getaccuracies.glmnet",
			"getaccuracies.pam",
			"getaccuracies.myDlda",
			"subgroupAcc",
			"sample_clusters")

	n_cores = future::availableCores()
	cl = makeCluster(n_cores)
	a=clusterCall(cl, function() library(glmnet))
	a=clusterCall(cl, function() library(pamr))
	a=clusterCall(cl, function() library(crc))
	a=clusterCall(cl, function() library(magrittr))
	clusterSetRNGStream(cl, iseed=123)
	clusterExport(cl=cl, varlist=wrappers)

	repl = parLapply(cl, 1:Nsim,
	function(i_sim, predictors, classes, id, GETACCURACIES, fctr, downsample, train_size) {
	    p = ncol(predictors)
	    samplep = sample(ncol(predictors),floor(ncol(predictors)*downsample))
 
	    	id_clusters = lapply(unique(id), function(x) which(id==x)) # id:{samples corresp to this id}
		id_class = sapply(unique(id), function(x) classes[which(id==x)][1]) # id:cluster
	    sc = sample_clusters(id_clusters, id_class, train_size)

	    	Ztrn = predictors[sc$samples_trn,samplep,drop=FALSE]
		Ztst = predictors[sc$samples_tst,samplep,drop=FALSE]
		Ttrn = classes[sc$samples_trn]
		Ttst = classes[sc$samples_tst]
		fctr_tst = fctr[sc$samples_tst]
		return(GETACCURACIES(Ztrn, Ttrn, Ztst, Ttst, fctr_tst))
	}, predictors=predictors, classes=classes, id=id, GETACCURACIES=GETACCURACIES, fctr=fctr,
	downsample=downsample, train_size=train_size)
	stopCluster(cl)
	
	n1=sum(classes==levels(classes)[1])
	n2=sum(classes==levels(classes)[2])
	n0=min(n1,n2)
	if (train_size < 1) n0trn = floor(n0*train_size)
	if (train_size >= 1) n0trn = train_size
    n0tst = n0 - n0trn
    p = ncol(predictors)
  
    results = list()
    results$out = sapply(repl, function(x) x$results) %>% t
    results$selected.ix = lapply(repl, function(x) x$selected.ix )
  
    results$clf = GETACCURACIES
    results$class_names = levels(classes)
    results$class_sizes = c(n1, n2)
    results$Nsim = Nsim
    results$n0 = n0
    results$n0tst = n0tst
    results$p = p
    save(results, file=filename)
}
