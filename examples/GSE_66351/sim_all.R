set.seed(123)
source("../../sim_params.R")
source("sim_fns.R")
clfs = c("crc", "glmnet", "pam", "myDlda")

# All cell types
data = readRDS("../../processed_data/GSE_66351/GSE66351_processed.RDS")
df = data$phenotype
predictors = data$beta
class = factor(df$`diagnosis:ch1`)
class1 = levels(class)[1]
class2 = levels(class)[2]
id = df$`donor_id:ch1` 			  # do not split clusters occuring across this factor
ct = factor(df$`cell type:ch1`) # report acc rates across this factor
br = factor(df$`brain_region:ch1`)
fctr = sapply(1:length(ct), function(i) ifelse(ct[i] == "bulk", paste(ct[i], br[i]), paste(ct[i]) ) )
fctr = factor(fctr)

for (clf in clfs) {
	filename = paste0(clf, '_',class1, "_vs_", class2, "_", ".out")
	run_sim_clustered(clf, predictors, class, id, fctr, Nsim, downsample = 1.0, 					train_size = 0.8, filename=filename)
}
