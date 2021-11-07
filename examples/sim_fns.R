crc_predict = function(fit, Ztst, clf = c("S", "L", "C")) {
  beta0 = fit[[paste0("beta0.", clf)]]
  beta1 = fit[[paste0("beta1.", clf)]]
  scores = Ztst %*% matrix(beta1) + beta0
  pred = ifelse(scores > 0, fit$class_names[2], fit$class_names[1])
  return(pred)
}

train_test_crc = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  min.eig.ratio = params$min.eig.ratio
  fit = crc(Ztrn, Ttrn, prior = c(0.5, 0.5), min.eig.ratio = min.eig.ratio)
  accS = mean(crc_predict(fit, Ztst, "S") == Ttst)
  accL = mean(crc_predict(fit, Ztst, "L") == Ttst)
  pred = crc_predict(fit, Ztst, "C")
  accC = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  # beta0 = fit[["beta0.C"]]
  # beta1 = fit[["beta1.C"]]
  # selected.ix = which(fit$S.fit$beta1.S.S != 0)
  # n.selected = length(selected.ix)
  return(list(accC = accC, accS = accS, accL = accL, subgroup_acc = subgroup_acc))
}

train_test_glmnet = function(Ztrn, Ttrn, Ztst, Ttst, params, fctr_tst=NULL) {
  alpha = params$alpha
  cv.fit = cv.glmnet(Ztrn, Ttrn, family="binomial", alpha=alpha, standardize=TRUE)
  pred = predict(cv.fit, newx=Ztst, type="class")
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  # beta1 = coef(cv.fit, s = "lambda.1se")[-1]
  # beta0 = coef(cv.fit, s = "lambda.1se")[1]
  # lambda = cv.fit$lambda.1se
  # bestlam = which(cv.fit$lambda==cv.fit$lambda.1se)
  # selected.ix = which(cv.fit$glmnet.fit$beta[,bestlam]!=0)
  # n.selected = length(selected.ix)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

train_test_pam = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  cv.fit = pamr.train(list("x"=t(Ztrn), "y"=Ttrn))
  t0 = max(which(cv.fit$errors == min(cv.fit$errors)))
  best.thresh = cv.fit$threshold[t0]
  selected.ix = pamr.listgenes(fit=cv.fit, list("x"=t(Ztrn), "y"=Ttrn, "geneid"=1:ncol(Ztrn)), threshold=best.thresh)[,"id"]
  pred = pamr.predict(cv.fit, newx=t(Ztst), threshold=best.thresh)
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  # n.selected = cv.fit$nonzero[t0]
  return(list(acc=acc, subgroup_acc = subgroup_acc))
}

myDlda = function(S, class, prior = c(0.5, 0.5), S.grid = crc:::get_S.grid("default", S)) {
  if (sum(apply(S, 2, sd) == 0) > 0) stop("One or more feature variances are constant.")
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

train_test_myDlda = function(Ztrn, Ttrn, Ztst, Ttst, params=NULL, fctr_tst=NULL) {
  fit = myDlda(Ztrn, Ttrn, prior = c(0.5, 0.5))
  beta1 = fit$beta1
  beta0 = fit$beta0
  scores = Ztst[, fit$selected.idx, drop = F] %*% fit$beta1 + fit$beta0
  pred = ifelse(scores > 0, fit$class_names[2], fit$class_names[1])
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  # n.selected = length(fit$selected.idx)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}

chunked_parLapply = function(nsim, f, params, predictors, classes, cl) {
  ncores = length(cl)
  # if (ncores == 1) warning("Only using one core in parallel execution.")
  x = 1:nsim
  nchunks = nsim / ncores
  # if (nchunks < 1) stop("ncores must be <= nsim.")
  splits = lapply(1:nchunks, function(i) rep(i, ncores))
  if (nsim %% floor(nchunks) != 0) {
    splits[[(nchunks+1)]] = rep((floor(nchunks)+1), nsim %% floor(nchunks))
  }
  i.list = split(x, factor(unlist(splits)))
  result.list = list()
  for (i in seq_along(i.list)) {
    i.vec = i.list[[i]]
    result.list[i.vec] = parLapply(cl, x[i.vec], 
                                   function(i_sim, params, predictors, classes) f(params, predictors, classes), 
                                   params = params, predictors = predictors, classes = classes)
  }
  return(result.list)
}


## Experiment (tracks runtimes)
experiment_rt = function(params, predictors, classes) {
  clf = match.fun(paste0("train_test_", params$clf))
  xp = unlist(params$xp, recursive = FALSE)
  train_size = params$train_size
  downsample = params$downsample
  
  n1 = sum(classes==levels(classes)[1])
  n2 = sum(classes==levels(classes)[2])
  n0 = min(n1,n2)
  if (train_size < 1) n0trn = floor(n0*train_size)
  if (train_size >= 1) n0trn = train_size
  if (n0trn >= n0) {
    print("train_size too large")
    return()
  }
  n0tst = n0 - n0trn
  p = ncol(predictors)
  
  sample1 = sample(which(classes==levels(classes)[1]), n0)
  sample2 = sample(which(classes==levels(classes)[2]), n0)
  samplep = sample(ncol(predictors),floor(ncol(predictors)*downsample))
  trn = c(sample1[1:n0trn], sample2[1:n0trn])
  tst = c(sample1[(n0trn+1):n0], sample2[(n0trn+1):n0])
  
  Ttrn = classes[trn]
  Ztrn = predictors[trn,samplep,drop=FALSE]
  Ttst = classes[tst]
  Ztst = predictors[tst,samplep,drop=FALSE]
  
  ## Drop constant features
  is_constant = apply(Ztrn, 2, function(x) sd(x)==0)
  if (sum(is_constant)!=0) cat("Dropping constant columns.\n")
  Ztrn = Ztrn[,!is_constant]
  Ztst = Ztst[,!is_constant]
  
  runtime = system.time({res = clf(Ztrn, Ttrn, Ztst, Ttst, xp)})['elapsed']
  res$runtime = runtime
  return(res)
}

## For leave-subject-out simulations, sample subjects (not observations)
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

## For leave-subject-out simulations, calculates accuracy within subgroups
subgroup_acc = function(pred, Ttst, fctr_tst) {
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

## Leave-subject-out experiment (tracks runtimes)
experiment_lso_rt = function(params, predictors, classes) {
  clf = match.fun(paste0("train_test_", params$clf))
  xp = unlist(params$xp, recursive = FALSE)
  train_size = params$train_size
  downsample = params$downsample
  id = unlist(params$id, recursive = FALSE)
  fctr = unlist(params$fctr, recursive = FALSE)
  
  samplep = sample(ncol(predictors),floor(ncol(predictors)*downsample))
  id_clusters = lapply(unique(id), function(x) which(id==x)) # id:{samples corresp to this id}
  id_class = sapply(unique(id), function(x) classes[which(id==x)][1]) # id:cluster
  sc = sample_clusters(id_clusters, id_class, train_size)
  Ztrn = predictors[sc$samples_trn,samplep,drop=FALSE]
  Ztst = predictors[sc$samples_tst,samplep,drop=FALSE]
  Ttrn = classes[sc$samples_trn]
  Ttst = classes[sc$samples_tst]
  fctr_tst = fctr[sc$samples_tst]
  
  ## Drop constant features
  is_constant = apply(Ztrn, 2, function(x) sd(x)==0)
  if (sum(is_constant)!=0) cat("Dropping constant columns.\n")
  Ztrn = Ztrn[,!is_constant]
  Ztst = Ztst[,!is_constant]
  
  runtime = system.time({res = clf(Ztrn, Ttrn, Ztst, Ttst, xp, fctr_tst)})['elapsed']
  res$runtime = runtime
  return(res)
}

## Function to run experiment
experiment = function(params, predictors, classes) {
  clf = match.fun(paste0("train_test_", params$clf))
  xp = unlist(params$xp, recursive = FALSE)
  train_size = params$train_size
  downsample = params$downsample

  n1 = sum(classes==levels(classes)[1])
  n2 = sum(classes==levels(classes)[2])
  n0 = min(n1,n2)
  if (train_size < 1) n0trn = floor(n0*train_size)
  if (train_size >= 1) n0trn = train_size
  if (n0trn >= n0) {
    print("train_size too large")
    return()
  }
  n0tst = n0 - n0trn
  p = ncol(predictors)

  sample1 = sample(which(classes==levels(classes)[1]), n0)
  sample2 = sample(which(classes==levels(classes)[2]), n0)
  samplep = sample(ncol(predictors),floor(ncol(predictors)*downsample))
  trn = c(sample1[1:n0trn], sample2[1:n0trn])
  tst = c(sample1[(n0trn+1):n0], sample2[(n0trn+1):n0])

  Ttrn = classes[trn]
  Ztrn = predictors[trn,samplep,drop=FALSE]
  Ttst = classes[tst]
  Ztst = predictors[tst,samplep,drop=FALSE]
  
  ## Drop constant features
  is_constant = apply(Ztrn, 2, function(x) sd(x)==0)
  if (sum(is_constant)!=0) cat("Dropping constant columns.\n")
  Ztrn = Ztrn[,!is_constant]
  Ztst = Ztst[,!is_constant]
  return(clf(Ztrn, Ttrn, Ztst, Ttst, xp))
}


## Function to run [l]eave [s]ubject [o]ut experiment
experiment_lso = function(params, predictors, classes) {
  clf = match.fun(paste0("train_test_", params$clf))
  xp = unlist(params$xp, recursive = FALSE)
  train_size = params$train_size
  downsample = params$downsample
  id = unlist(params$id, recursive = FALSE)
  fctr = unlist(params$fctr, recursive = FALSE)

  samplep = sample(ncol(predictors),floor(ncol(predictors)*downsample))
  id_clusters = lapply(unique(id), function(x) which(id==x)) # id:{samples corresp to this id}
  id_class = sapply(unique(id), function(x) classes[which(id==x)][1]) # id:cluster
  sc = sample_clusters(id_clusters, id_class, train_size)
  Ztrn = predictors[sc$samples_trn,samplep,drop=FALSE]
  Ztst = predictors[sc$samples_tst,samplep,drop=FALSE]
  Ttrn = classes[sc$samples_trn]
  Ttst = classes[sc$samples_tst]
  fctr_tst = fctr[sc$samples_tst]
  
  ## Drop constant features
  is_constant = apply(Ztrn, 2, function(x) sd(x)==0)
  if (sum(is_constant)!=0) cat("Dropping constant columns.\n")
  Ztrn = Ztrn[,!is_constant]
  Ztst = Ztst[,!is_constant]
  
  return(clf(Ztrn, Ttrn, Ztst, Ttst, xp, fctr_tst))
}