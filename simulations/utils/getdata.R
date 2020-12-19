getdata = function(n,p) {
  nb = 3
  k = 3
  T = matrix(rep(c(-1, 1), each=n/2))
  L1 = matrix(rnorm(n*k),n,k)
  alpha = matrix(rnorm(p*k),k,p)
  beta = matrix(0,1,p); beta[1:nb] = 1/sqrt(nb)
  epsilon = matrix(rnorm(n*p),n,p)
  S = T%*%beta + epsilon
  Z1 = S + L1%*%alpha
  eta = matrix(1,1,k)/sqrt(k)
  L2 = L1 + T%*%eta
  Z2 = S + L2%*%alpha
  return(list(T=T,S=S,Z1=Z1,Z2=Z2,alpha=alpha,beta=beta,eta=eta))
}

