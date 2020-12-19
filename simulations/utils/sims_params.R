# Simulation parameters

param = data.frame(rbind(
cbind(
  Nsim = c(100, 100, 100, 30, 10),
  n = c(50, 100, 200, 500, 1000),
  p = 20000),
cbind(
  Nsim = c(100, 100, 100, 30, 10),
  n = c(50, 100, 200, 500, 1000),
  p = 100000),
cbind(
  Nsim = c(100, 100, 100, 30, 10),
  n = c(50, 100, 200, 500, 1000),
  p = 500000)
))

param_p20k = data.frame(cbind(
  Nsim = c(100, 100, 100, 30, 10),
  n = c(50, 100, 200, 500, 1000),
  p = 20000))

param_p100k = data.frame(cbind(
		Nsim = c(100, 100, 100, 30, 10),
		n = c(50, 100, 200, 500, 1000),
		p = 100000))
		
param_p500k = data.frame(cbind(
		Nsim = c(100, 100, 100, 30, 10),
		n = c(50, 100, 200, 500, 1000),
		p = 500000
))