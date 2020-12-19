set.seed(123)
source("../sim_params.R")
source("../sim_fns.R")
clfs = c("crc", "glmnet", "pam", "myDlda")

# Data
data = readRDS("../processed_data/E_MTAB_1532/EMTAB1532_processed.RDS")
predictors = data$Y

# Simulation
classes = data$phenotype$Characteristics.disease.state.
class1 = levels(factor(classes))[1]
class2 = levels(factor(classes))[2]
for (clf in clfs) {
    filename = paste0(clf, '_',class1, "_vs_", class2, ".out")
	run_sim(clf, predictors, classes, Nsim, downsample = 1.0, train_size = 0.8, filename)
}