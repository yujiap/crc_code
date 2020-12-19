library(glmnet)
library(magrittr)

getaccuracies = function(sim, alpha) {
  car = function(a, b, m, v) {
    mm1 = a - m%*%b
    mp1 = a + m%*%b
    #v = t(b)%*%sigma%*%b
    accm1 = pnorm(0, mean=mm1, sd=sqrt(v), lower.tail=TRUE)
    accp1 = pnorm(0, mean=mp1, sd=sqrt(v), lower.tail=FALSE)
    return(0.5*(accm1+accp1))
  }
  results = matrix(NA,3,3)
  colnames(results) = c("glmnet", "lambda.1se", "Nsel")
  rownames(results) = c("S", "L1", "L2")
  
  # Simple
  m = sim$beta
  cv.fit = cv.glmnet(sim$S, sim$T, family="binomial", alpha=alpha, standardize=T)
  beta1 = coef(cv.fit, s = "lambda.1se")[-1]
  beta0 = coef(cv.fit, s = "lambda.1se")[1]
  v = t(beta1) %*% beta1
  results[1,1] = car(beta0, beta1, m, v)
  results[1,2] = cv.fit$lambda.1se
  results[1,3] = sum(beta1!=0)
  
  # Uncorrelated
  m = sim$beta
  cv.fit = cv.glmnet(sim$Z1, sim$T, family="binomial", alpha=alpha, standardize=T)
  beta1 = coef(cv.fit, s = "lambda.1se")[-1]
  beta0 = coef(cv.fit, s = "lambda.1se")[1]
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[2,1] = car(beta0, beta1, m, v)
  results[2,2] = cv.fit$lambda.1se
  results[2,3] = sum(beta1!=0)
  
  # Correlated
  m = sim$beta + sim$eta%*%sim$alpha
  cv.fit = cv.glmnet(sim$Z2, sim$T, family="binomial", alpha=alpha, standardize=T)
  beta1 = coef(cv.fit, s = "lambda.1se")[-1]
  beta0 = coef(cv.fit, s = "lambda.1se")[1]
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[3,1] = car(beta0, beta1, m, v)
  results[3,2] = cv.fit$lambda.1se
  results[3,3] = sum(beta1!=0)
  return(results)    
}


print_summary = function(results) {
  n = results$n
  p = results$p
  accuracies = results$accuracies
  Nsim = dim(accuracies)[3]
  cat("\n\nAlpha: ", alpha, "\n")
  cat("\n\nDimensions: ", c(n,p), "\n")
    cat("\nAvg results: \n")
    out = cbind(apply(accuracies,c(1,2),mean)[,1:2],
                apply(accuracies, 3, function(mtx) mtx[,3]) %>% apply(., 1, median))
    colnames(out)[3] = "n.selected"
    print(out)
    cat("SE: \n")
    out = apply(accuracies,c(1,2),sd)/sqrt(Nsim)
    out[,3] = out[,3] * 1.25
    print(out)
}