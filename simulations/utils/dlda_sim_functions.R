library(crc)
library(magrittr)

# DLDA with feature selection
myDlda = function(S, class, prior = c(0.5, 0.5), S.grid = crc:::get_S.grid("default", S)) {
  n = length(class)
  n1 = sum(class == -1)
  n2 = sum(class == 1)
  I = crc:::get_I(class)
  Cs = t(I) %*% S # [C]ol[s]ums of S
  Cs.sq = t(I) %*% S^2 # [C]ol[s]ums of S [sq]uared
  Mu1 = Cs[1,] / n1
  Mu2 = Cs[2,] / n2
  ss1 = Cs.sq[1,] - Cs[1,]^2/n1
  ss2 = Cs.sq[2,] - Cs[2,]^2/n2
  
  nN = length(S.grid)
  scores.LOO_N = matrix(0, n, nN)
  for (i in 1:n) { # check which is faster, for vs. sapply. also beware scope
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
                    sm1 = mean(x[class == -1])
                    sm2 = mean(x[class == 1])
                    spv = (n1*var(x[class == -1]) + n2*var(x[class == 1]))/(n-2)
                    p.mistake = prior[1]*pnorm(sm1/sqrt(spv)) + prior[2]*pnorm(-sm2/sqrt(spv))
                    acc = 1 - p.mistake
                  })
  
  sigmas2 = (ss1 + ss2) / (n-2)
  t = (Mu2 - Mu1) / sqrt(sigmas2 * n/(n1*n2)) # recheck
  
  idx = which.max(thr_acc)
  topN = order(abs(t), decreasing = T)[1:S.grid[idx]]
  beta1 = (Mu2[topN] - Mu1[topN]) / sigmas2[topN]
  beta0 = log(prior[2]/prior[1]) - 0.5*sum((Mu2[topN] + Mu1[topN])*beta1)
  
  return(list("beta1" = beta1,
              "beta0" = beta0,
              "scores.LOO"=scores.LOO_N[,idx],
              "selected.idx"=topN,
              "v"=thr_acc
  ))
}


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
  colnames(results) = c("dlda", "best.thresh", "Nsel")
  rownames(results) = c("S", "L1", "L2")
  
  # Simple
  m = sim$beta
  p = ncol(sim$S)
  fit = myDlda(sim$S, sim$T)
  beta1 = rep(0, p); beta1[fit$selected.idx] = fit$beta1
  beta0 = fit$beta0
  v = t(beta1) %*% beta1
  results[1,1] = car(beta0, beta1, m, v)
  results[1,2] = 0
  results[1,3] = sum(beta1!=0)
  
  # Uncorrelated
  m = sim$beta
  fit = myDlda(sim$Z1, sim$T)
  beta1 = rep(0, p); beta1[fit$selected.idx] = fit$beta1
  beta0 = fit$beta0
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[2,1] = car(beta0, beta1, m, v)
  results[2,2] = 0
  results[2,3] = sum(beta1!=0)
  
  # Correlated
  m = sim$beta + sim$eta%*%sim$alpha
  fit = myDlda(sim$Z2, sim$T)
  beta1 = rep(0, p); beta1[fit$selected.idx] = fit$beta1
  beta0 = fit$beta0
  v = t(beta1) %*% beta1 + sum((t(beta1) %*% t(sim$alpha))^2)
  results[3,1] = car(beta0, beta1, m, v)
  results[3,2] = 0
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