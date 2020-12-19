library(crc)
library(magrittr)

getaccuracies = function(sim) {
  car = function(a, b, m, v) {
    mm1 = a - m%*%b
    mp1 = a + m%*%b
    #v = t(b)%*%sigma%*%b
    accm1 = pnorm(0, mean=mm1, sd=sqrt(v), lower.tail=TRUE)
    accp1 = pnorm(0, mean=mp1, sd=sqrt(v), lower.tail=FALSE)
    return(0.5*(accm1+accp1))
  }
  p = ncol(sim$S)
  results = matrix(NA,3,4)
  colnames(results) = c("S", "L", "C", "Nsel")
  rownames(results) = c("S", "L1", "L2")
  
  
  ## Simple
  m = sim$beta
  # sigma = diag(p)
  crc = crc(Z = sim$S, classes = sim$T, min.eig.ratio = 0)
  v = t(crc$beta1.S) %*% crc$beta1.S
  results[1,1] = car(crc$beta0.S, crc$beta1.S, m, v)
  v = t(crc$beta1.L) %*% crc$beta1.L
  results[1,2] = car(crc$beta0.L, crc$beta1.L, m, v)
  v = t(crc$beta1.C) %*% crc$beta1.C
  results[1,3] = car(crc$beta0.C, crc$beta1.C, m, v)
  results[1,4] = length(crc$S.selected)
  
  
  ## Uncorrelated
  m = sim$beta
  # sigma = diag(p) + t(sim$alpha)%*%sim$alpha
  crc = crc(Z = sim$Z1, classes = sim$T, min.eig.ratio = 0)
  v = t(crc$beta1.S) %*% crc$beta1.S + sum((t(crc$beta1.S) %*% t(sim$alpha))^2)
  results[2,1] = car(crc$beta0.S, crc$beta1.S, m, v)
  v = t(crc$beta1.L) %*% crc$beta1.L + sum((t(crc$beta1.L) %*% t(sim$alpha))^2)    
  results[2,2] = car(crc$beta0.L, crc$beta1.L, m, v)
  v = t(crc$beta1.C) %*% crc$beta1.C + sum((t(crc$beta1.C) %*% t(sim$alpha))^2)    
  results[2,3] = car(crc$beta0.C, crc$beta1.C, m, v)
  results[2,4] = length(crc$S.selected)
  
  
  # Correlated
  m = sim$beta + sim$eta%*%sim$alpha
  #sigma = diag(p) + t(sim$alpha)%*%sim$alpha
  crc = crc(Z = sim$Z2, classes = sim$T, min.eig.ratio = 0)
  v = t(crc$beta1.S) %*% crc$beta1.S + sum((t(crc$beta1.S) %*% t(sim$alpha))^2)
  results[3,1] = car(crc$beta0.S, crc$beta1.S, m, v)
  v = t(crc$beta1.L) %*% crc$beta1.L + sum((t(crc$beta1.L) %*% t(sim$alpha))^2)    
  results[3,2] = car(crc$beta0.L, crc$beta1.L, m, v)
  v = t(crc$beta1.C) %*% crc$beta1.C + sum((t(crc$beta1.C) %*% t(sim$alpha))^2)    
  results[3,3] = car(crc$beta0.C, crc$beta1.C, m, v)
  results[3,4] = length(crc$S.selected)
  
  return(results)    
}


print_summary = function(results) {
  n = results$n
  p = results$p
  accuracies = results$accuracies
  Nsim = dim(accuracies)[3]

  cat("\n\nDimensions: ", c(n,p), "\n")
  cat("\nAvg results: \n")
  out = cbind(apply(accuracies,c(1,2),mean)[,1:3],
              apply(accuracies, 3, function(mtx) mtx[,4]) %>% apply(., 1, median))
  colnames(out)[4] = "Nsel"
  print(out)
  
  cat("SE: \n")
  out = apply(accuracies, c(1,2),sd)/sqrt(Nsim)
  out[,3] = out[,3] * 1.25
  print(out)
}