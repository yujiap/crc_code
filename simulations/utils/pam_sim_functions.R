library(pamr)
library(magrittr)

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
    xbar.k = cv.fit$centroids[,lab]
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

# -----------------------------------

getaccuracies = function(sim) {
  car = function(a, b, m, v) {
    mm1 = a - m%*%b
    mp1 = a + m%*%b
    #v = t(b)%*%sigma%*%b
    accm1 = pnorm(0, mean=mm1, sd=sqrt(v), lower.tail=TRUE)
    accp1 = pnorm(0, mean=mp1, sd=sqrt(v), lower.tail=FALSE)
    return(0.5*(accm1+accp1))
  }
  results = matrix(NA,3,3)
  colnames(results) = c("pam", "best.thresh", "Nsel")
  rownames(results) = c("S", "L1", "L2")
  
  # Simple
  m = sim$beta
  cv.fit = pamr.train(list("x"=t(sim$S), "y"=sim$T))
  pam_est = pamEstimates(cv.fit)
  pam_coef = pamCoef(pam_est, prior = c(0.5, 0.5))
  beta1 = pam_coef$beta1
  beta0 = pam_coef$beta0
  v = t(beta1) %*% beta1
  results[1,1] = car(beta0, beta1, m, v)
  results[1,2] = pam_est$best.thresh
  results[1,3] = sum(beta1!=0)
  
  # Uncorrelated
  m = sim$beta
  cv.fit = pamr.train(list("x"=t(sim$Z1), "y"=sim$T))
  pam_est = pamEstimates(cv.fit)
  pam_coef = pamCoef(pam_est, prior = c(0.5, 0.5))
  beta1 = pam_coef$beta1
  beta0 = pam_coef$beta0
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[2,1] = car(beta0, beta1, m, v)
  results[2,2] = pam_est$best.thresh
  results[2,3] = sum(beta1!=0)
  
  # Correlated
  m = sim$beta + sim$eta%*%sim$alpha
  cv.fit = pamr.train(list("x"=t(sim$Z2), "y"=sim$T))
  pam_est = pamEstimates(cv.fit)
  pam_coef = pamCoef(pam_est, prior = c(0.5, 0.5))
  beta1 = pam_coef$beta1
  beta0 = pam_coef$beta0
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[3,1] = car(beta0, beta1, m, v)
  results[3,2] = pam_est$best.thresh
  results[3,3] = sum(beta1!=0)
  return(results)    
}


print_summary = function(results) {
  n = results$n
  p = results$p
  accuracies = results$accuracies
  Nsim = dim(accuracies)[3]

  cat("\n\nDimensions: ", c(n,p), "\n")
  cat("\nAvg results: \n")
  out = cbind(apply(accuracies,c(1,2),mean)[,1:2],
              apply(accuracies, 3, function(mtx) mtx[,3]) %>% apply(., 1, median))
  colnames(out)[3] = "Nsel"
  print(out)
  cat("SE: \n")
  out = apply(accuracies,c(1,2),sd)/sqrt(Nsim)
  out[,3] = out[,3] * 1.25
  print(out)
}