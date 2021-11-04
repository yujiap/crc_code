# sample_gaussian_mixture = function(n) {
#   r = 3
#   p1 = 0.2
#   p2 = 0.3
#   p3 = 0.5
#   mu1 = c(0, 0, 0)
#   mu2 = c(0.2, 0.2, 0.2)
#   mu3 = (p1*mu1 + p2*mu2) / -p3
#   c1 = diag(r)
#   c2 = diag(r)
#   c3 = (p1*c1 + p2*c2 + p1*tcrossprod(mu1) + p2*tcrossprod(mu2) + p3*tcrossprod(mu3) - diag(r)) / (-p3)
#   id = sample(1:3, prob=c(p1, p2, p3), size=n, replace=TRUE)
#   mus = list(mu1, mu2, mu3)
#   covs = list(c1, c2, c3)
#   samples = t(sapply(id, function(x) mvrnorm(1, mu = mus[x][[1]], Sigma = covs[x][[1]])))
#   return(samples)
# }

sample_univariate_gaussian_mixture = function(n) {
  p1 = 0.2
  p2 = 0.3
  p3 = 0.5
  mu1 = 1.5
  mu2 = -0.5
  mu3 = (p1*mu1 + p2*mu2) / -p3
  c1 = 0.02
  c2 = 0.6
  c3 = (p1*c1 + p2*c2 + p1*as.numeric(tcrossprod(mu1)) + p2*as.numeric(tcrossprod(mu2)) + p3*as.numeric(tcrossprod(mu3)) - 1) / (-p3)
  id = sample(1:3, prob=c(p1, p2, p3), size=n, replace=TRUE)
  mus = list(mu1, mu2, mu3)
  covs = list(c1, c2, c3)
  print(mus)
  print(covs)
  samples = as.vector(sapply(id, function(x) mvrnorm(1, mu = mus[x][[1]], Sigma = covs[x][[1]])))
  return(samples)
}

get_data = function(n, p, frac1 = 0.5, alpha.rs.level = 1.0, L.type = "normal", eps.type = "normal") {
  nb = 3
  r = 3
  T = matrix(rep(1,n))
  T[1:floor(n*frac1)] = -1
  alpha = t(sapply(1:r, function(l) {
    DE = sample(1:p, floor(alpha.rs.level*p))
    alpha.i = rep(0, p)
    alpha.i[DE] = rnorm(length(DE))
    alpha.i
  }))
  gamma = matrix(0, 1, p); gamma[1:nb] = 1/sqrt(nb)
  eta = matrix(1, 1, r)/sqrt(r)
  if (eps.type == "normal") {
    epsilon = matrix(rnorm(n*p), n, p)
  } else if (eps.type == "dexp") {
    sc = 4
    epsilon = matrix(rdexp(n*p, location = 0, scale = sc), n, p)
    epsilon = epsilon / sqrt(2 * sc^2)
  } else if (eps.type == "unif") {
    epsilon = matrix(runif(n*p, min = 0, max = 1), n, p)
    epsilon = (epsilon - 0.5) * sqrt(12)
  } else if (eps.type == "t") {
    ndf = 3 # > 2
    epsilon = matrix(rt(n*p, ndf), n, p) / sqrt(ndf/(ndf-2))
  }
  if (L.type == "normal") {
    L1 = matrix(rnorm(n*r), n, r)
  } else if (L.type == "mixture") {
    L1 = sapply(1:r, function(ir) sample_univariate_gaussian_mixture(n))
  } else if (L.type == "uniform") {
    L1 = matrix(runif(n*r, min = 0, max = 1), nrow =  n, ncol = r)
    L1 = (L1 - 0.5) * sqrt(12)
  } else if (L.type == "skewnormal") {
    sn_alpha = 10
    sn_delta = sn_alpha/sqrt(1+sn_alpha^2)
    sn_omega = sqrt(pi/(pi-2*sn_delta^2))
    sn_xi = -sn_omega*sn_delta*sqrt(2/pi)
    L1 = matrix(rsn(n=n*r, xi=sn_xi, omega=sn_omega, alpha=sn_alpha, tau=0), n, r)
  }
  S = T%*%gamma + epsilon
  Z1 = S + L1%*%alpha
  L2 = L1 + T%*%eta
  Z2 = S + L2%*%alpha
  return(list(T=T,S=S,Z1=Z1,Z2=Z2,alpha=alpha,gamma=gamma,eta=eta))
}

generate_epsilon = function(n, p, eps.type = "normal") {
  if (eps.type == "normal") {
    epsilon = matrix(rnorm(n*p), n, p)
  } else if (eps.type == "dexp") {
    sc = 4
    epsilon = matrix(rdexp(n*p, location = 0, scale = sc), n, p)
    epsilon = epsilon / sqrt(2 * sc^2)
  } else if (eps.type == "unif") {
    epsilon = matrix(runif(n*p, min = 0, max = 1), n, p)
    epsilon = (epsilon - 0.5) * sqrt(12)
  } else if (eps.type == "t") {
    ndf = 3 # > 2
    epsilon = matrix(rt(n*p, ndf), n, p) / sqrt(ndf/(ndf-2))
  } else {
    stop("Invalid eps.type")
  }
  return(epsilon)
}

generate_L = function(n, r, L.type = "normal") {
  if (L.type == "normal") {
    L1 = matrix(rnorm(n*r), n, r)
  } else if (L.type == "mixture") {
    L1 = sapply(1:r, function(ir) sample_univariate_gaussian_mixture(n))
  } else if (L.type == "uniform") {
    L1 = matrix(runif(n*r, min = 0, max = 1), nrow =  n, ncol = r)
    L1 = (L1 - 0.5) * sqrt(12)
  } else if (L.type == "skewnormal") {
    sn_alpha = 10
    sn_delta = sn_alpha/sqrt(1+sn_alpha^2)
    sn_omega = sqrt(pi/(pi-2*sn_delta^2))
    sn_xi = -sn_omega*sn_delta*sqrt(2/pi)
    L1 = matrix(rsn(n=n*r, xi=sn_xi, omega=sn_omega, alpha=sn_alpha, tau=0), n, r)
  } else {
    stop("Invalid L.type")
  }
  return(L1)
}

generate_alpha = function(p, r, alpha.rs.level = 1.0) {
  alpha = t(sapply(1:r, function(l) {
    DE = sample(1:p, floor(alpha.rs.level*p))
    alpha.i = rep(0, p)
    alpha.i[DE] = rnorm(length(DE))
    alpha.i
  }))
  return(alpha)
}

generate_T = function(n, frac1) {
  T = matrix(rep(1, n))
  T[1:floor(n*frac1)] = -1
  return(T)
}

get_data_trn_tst = function(n, ntst, p, frac1 = 0.5, alpha.rs.level = 1.0, L.type = "normal", eps.type = "normal") {
  nb = 3
  r = 3
  T = generate_T(n, frac1)
  Ttst = generate_T(ntst, frac1)
  alpha = generate_alpha(p, r, alpha.rs.level = alpha.rs.level)
  gamma = matrix(0, 1, p); gamma[1:nb] = 1/sqrt(nb)
  eta = matrix(1, 1, r)/sqrt(r)
  epsilon = generate_epsilon(n, p, eps.type = eps.type)
  epsilon_tst = generate_epsilon(ntst, p, eps.type = eps.type)
  L1 = generate_L(n, r, L.type = L.type)
  L1tst = generate_L(ntst, r, L.type = L.type)
  
  S = T%*%gamma + epsilon
  Z1 = S + L1%*%alpha
  L2 = L1 + T%*%eta
  Z2 = S + L2%*%alpha
  
  Stst = Ttst%*%gamma + epsilon_tst
  Z1tst = Stst + L1tst%*%alpha
  L2tst = L1tst + Ttst%*%eta
  Z2tst = Stst + L2tst%*%alpha
  return(list(T=T, S=S, Z1=Z1, Z2=Z2, 
              Ttst=Ttst, Stst=Stst, Z1tst=Z1tst, Z2tst=Z2tst, 
              alpha=alpha, gamma=gamma, eta=eta))
}

get_theoretical_accuracy = function(data, model, beta1, beta0) {
  car = function(a, b, m, v) {
    mm1 = a - m%*%b
    mp1 = a + m%*%b
    # v = t(b) %*% sigma %*% b
    accm1 = pnorm(0, mean=mm1, sd=sqrt(v), lower.tail=TRUE)
    accp1 = pnorm(0, mean=mp1, sd=sqrt(v), lower.tail=FALSE)
    return(0.5*(accm1+accp1))
  }
  if (model == "simple") {
    m = data$gamma
    v = t(beta1) %*% beta1
    # sigma = diag(p)
  } else if (model == "uncorrelated") {
    m = data$gamma
    v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(data$alpha))^2)
    # sigma = diag(p) + t(sim$alpha)%*%sim$alpha
  } else if (model == "correlated") {
    m = data$gamma + data$eta %*% data$alpha
    v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(data$alpha))^2)
    # sigma = diag(p) + t(sim$alpha)%*%sim$alpha
  }
  return(car(beta0, beta1, m, v))
}

chunked_parLapply = function(nsim, f, params, ncores, cl) {
  x = 1:nsim
  nchunks = nsim / ncores
  splits = lapply(1:nchunks, function(i) rep(i, ncores))
  if (nsim %% floor(nchunks) != 0) {
    splits[[(nchunks+1)]] = rep((floor(nchunks)+1), nsim %% floor(nchunks))
  }
  i.list = split(x, factor(unlist(splits)))
  result.list = list()
  for (i in seq_along(i.list)) {
    i.vec = i.list[[i]]
    result.list[i.vec] = parLapply(cl, x[i.vec], function(i_sim, params) f(params), params = params)
  }
  return(result.list)
}

## Sim functions for PAM

posPart = function(x) {
  ret = ifelse(x > 0, x, 0)
  return(ret)
}

# Select greatest threshold giving min error
pamEstimates = function(cv.fit) {
  best.thresh = which(cv.fit$errors == min(cv.fit$errors)) %>% max %>% cv.fit$threshold[.]
  ret = list()
  ret$var = (cv.fit$sd)^2
  for (lab in c("-1", "1")) {
    xbar = cv.fit$centroid.overall
    xbar.k = cv.fit$centroids[, lab]
    m.k = cv.fit$se.scale[lab]
    d.k = (xbar.k - xbar) / (m.k * cv.fit$sd)
    d.k_shrink = sign(d.k) * posPart(abs(d.k) - best.thresh)
    ret[[paste0("mean", lab)]] = d.k_shrink
  }
  ret$best.thresh = best.thresh
  return(ret)
}

pamCoef = function(pam_est, prior) {
  ret = list()
  ret$beta1 = (pam_est$mean1 - pam_est$`mean-1`) / pam_est$var
  ret$beta0 = -0.5 * sum((pam_est$mean1 + pam_est$`mean-1`) * ret$beta1) + log(prior[2]/prior[1])
  return(ret)
}

## DLDA with feature selection
myDlda <- function(S, classes, prior = c(0.5, 0.5), S.grid = crc:::get_S.grid("default", S)) {
  n <- length(classes)
  I <- crc:::get_I(classes)
  n1 <- sum(I[,1])
  n2 <- sum(I[,2])
  Cs <- t(I) %*% S # [C]ol[s]ums of S
  Cs.sq <- t(I) %*% S^2 # [C]ol[s]ums of S [sq]uared
  Mu1 <- Cs[1,] / n1
  Mu2 <- Cs[2,] / n2
  ss1 <- Cs.sq[1,] - Cs[1,]^2/n1
  ss2 <- Cs.sq[2,] - Cs[2,]^2/n2
  
  nN <- length(S.grid)
  scores.LOO_N <- matrix(0, n, nN)
  for (i in 1:n) { # check which is faster, for vs. sapply. also beware scope
    sigmas2 <- crc:::downdatePooledVar(S, I, Cs, Cs.sq, ss1, ss2, i)
    mu1 <- if (I[i, 1]) (Cs[1,] - S[i,])/(n1-1) else Mu1
    mu2 <- if (I[i, 2]) (Cs[2,] - S[i,])/(n2-1) else Mu2
    n1.n2 <- colSums(I[-i,])
    tdrop <- ((mu2 - mu1) / sqrt(sigmas2)) * sqrt((n1.n2[1]*n1.n2[2])/(n-1))
    toporder <- order(abs(tdrop), decreasing = TRUE)
    for (j in 1:nN) {
      N <- S.grid[j]
      topN <- toporder[1:N]
      beta1 <- (mu2[topN] - mu1[topN]) / sigmas2[topN]
      beta0 <- log(prior[2]/prior[1]) - 0.5 * sum((mu2[topN] + mu1[topN])*beta1)
      scores.LOO_N[i,j] <- S[i, topN, drop = F] %*% beta1 + beta0
    }
  }
  
  thr_acc <- apply(scores.LOO_N, 2,
                   function(x) {
                     sm1 <- mean(x[classes == -1])
                     sm2 <- mean(x[classes == 1])
                     spv <- (n1*var(x[classes == -1]) + n2*var(x[classes == 1]))/(n-2)
                     p.mistake <- prior[1]*pnorm(sm1/sqrt(spv)) + prior[2]*pnorm(-sm2/sqrt(spv))
                     acc <- 1 - p.mistake
                   })
  
  sigmas2 <- (ss1 + ss2) / (n-2)
  t <- (Mu2 - Mu1) / sqrt(sigmas2 * n/(n1*n2)) # recheck
  
  idx <- which.max(thr_acc)
  topN <- order(abs(t), decreasing = TRUE)[1:S.grid[idx]]
  beta1 <- rep(0, ncol(S))
  beta1[topN] <- (Mu2[topN] - Mu1[topN]) / sigmas2[topN]
  beta0 <- log(prior[2]/prior[1]) - 0.5*sum((Mu2[topN] + Mu1[topN])*beta1[topN])
  
  return(list("beta1" = beta1,
              "beta0" = beta0,
              "scores.LOO"=scores.LOO_N[,idx],
              "selected.idx"=topN,
              "v"=thr_acc
  ))
}

## [For simulations] Wrapper functions for non-CRC classifiers

wrapper_glmnet = function(Ztrn, Ttrn, alpha = 1) {
  # Use default 10 fold CV if # observations per fold is at least 5. Otherwise use LOOCV.
  nf = if (length(Ttrn)/10 < 5) length(Ttrn) else 10
  cv.fit = cv.glmnet(Ztrn, Ttrn, family="binomial", alpha=alpha, standardize=TRUE, nfolds = nf) ##
  beta1 = coef(cv.fit, s = "lambda.1se")[-1]
  beta0 = coef(cv.fit, s = "lambda.1se")[1]
  lambda = cv.fit$lambda.1se
  bestlam = which(cv.fit$lambda==cv.fit$lambda.1se)
  selected.ix = which(cv.fit$glmnet.fit$beta[,bestlam]!=0)
  n.selected = length(selected.ix)
  return(list(beta1=beta1, beta0=beta0, lambda=lambda, n.selected=n.selected, selected.ix=selected.ix))
}

wrapper_pam = function(Ztrn, Ttrn) {
  cv.fit = pamr.train(list("x"=t(Ztrn), "y"=Ttrn))
  pam_est = pamEstimates(cv.fit)
  pam_coef = pamCoef(pam_est, prior = c(0.5, 0.5))
  beta1 = pam_coef$beta1
  beta0 = pam_coef$beta0
  t0 = which(cv.fit$errors == min(cv.fit$errors)) %>% max
  selected.ix = which(pam_coef$beta1!=0)
  n.selected = cv.fit$nonzero[t0]
  best.thresh = cv.fit$threshold[t0]
  return(list(beta1=beta1, beta0=beta0, best.thresh=best.thresh, n.selected=n.selected, selected.ix=selected.ix))
}

wrapper_myDlda = function(Ztrn, Ttrn) {
  fit = myDlda(Ztrn, Ttrn, prior = c(0.5, 0.5))
  beta1 = fit$beta1
  beta0 = fit$beta0
  selected.ix = fit$selected.idx
  n.selected = length(selected.ix)
  return(list(beta1=beta1, beta0=beta0, n.selected=n.selected, selected.ix=selected.ix))
}

## Train/test functions

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

myDlda_trn_tst = function(S, class, prior = c(0.5, 0.5), S.grid = crc:::get_S.grid("default", S)) {
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
  fit = myDlda_trn_tst(Ztrn, Ttrn, prior = c(0.5, 0.5))
  beta1 = fit$beta1
  beta0 = fit$beta0
  scores = Ztst[, fit$selected.idx, drop = F] %*% fit$beta1 + fit$beta0
  pred = ifelse(scores > 0, fit$class_names[2], fit$class_names[1])
  acc = mean(pred == Ttst)
  subgroup_acc = if (!is.null(fctr_tst)) subgroup_acc(pred, Ttst, fctr_tst)
  # n.selected = length(fit$selected.idx)
  return(list(acc = acc, subgroup_acc = subgroup_acc))
}
