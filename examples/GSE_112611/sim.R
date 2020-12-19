set.seed(123)
source("../sim_params.R")
source("../sim_fns.R")
clfs = c("crc", "glmnet", "pam", "myDlda")

# Data
data = readRDS("../processed_data/GSE_112611/GSE112611_processed.RDS")
BL.ix = which(data$phenotype$`baseline vs follow-up:ch1` == "BL")
predictors = data$Y[BL.ix,]
class = factor(data$phenotype$`diagnosis:ch1`)[BL.ix]
tab = table(class)
nclass = length(tab)

# Simulation
system.time(
  for (i in 1:(nclass-1)) {
    for (j in (i+1):nclass) {
      class1 = names(tab[i])
      class2 = names(tab[j])
      ix = class %in% c(class1, class2)
      predictors_ix = subset(predictors, ix)
      classes_ix = class[ix]
      for (clf in clfs) {
        print(clf)
        filename = paste0(clf, '_',class1, "_vs_", class2, ".out")
        run_sim(clf, predictors_ix, classes_ix, Nsim, downsample = 1.0, train_size = 0.8, filename)
      }
    }
  }
)