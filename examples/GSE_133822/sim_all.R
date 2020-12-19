set.seed(123)
source("../../sim_params.R")
source("sim_fns.R")
clfs = c("crc", "glmnet", "pam", "myDlda")

data = readRDS("../../processed_data/GSE_133822/GSE133822_processed.RDS")
predictors = data$Y
df = data$phenotype

# Drop "Cancer" samples (very few per cell once stratified by cell type)
cancer_ix = which(df$patient_group=="Cancer")
predictors = predictors[-cancer_ix,]
df = df[-cancer_ix,]
class = factor(df$patient_group)
nclass = length(table(class))

for (c1 in 1:(nclass-1)) {
	for (c2 in (c1+1):nclass) {
		class1 = levels(class)[c1]
		class2 = levels(class)[c2]
		ix = (class==class1 | class==class2)
		for (clf in clfs) {
			filename = paste0(clf, '_',class1, "_vs_", class2, ".out")
			run_sim_clustered(clf, 
				predictors = predictors[ix,], 
				class = factor(df$patient_group[ix]), 
				id = df$subject_id[ix], # do not split clusters occuring across this factor
				fctr = factor(df$cell_type[ix]), # report acc rates across this factor
				Nsim, 
				downsample = 1.0, train_size = 0.8, filename=filename)
		}		
	}
}
