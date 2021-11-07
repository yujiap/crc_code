library(future)
library(parallel)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(here)
source(here("sim_fns.R"))
source(here("additional_methods", "sim_functions.R"))

## Setup parallel exec
n_cores = 10 # future::availableCores()
cl = makeCluster(n_cores)
packages = c("magrittr", "caret", "here")
a = clusterExport(cl, list("packages"))
a = clusterCall(cl, function() lapply(packages, library, character.only = TRUE))
a = clusterCall(cl, function() source(here("sim_fns.R")))
a = clusterCall(cl, function() source(here("additional_methods", "sim_functions.R")))
clusterSetRNGStream(cl, iseed = 123)

## Setup simulation (xp: e[x]tra [p]arameters)
sim = crossing(dataset = "GSE_66351",
               clf = "SVM",
               acc = NA,
               subgroup_acc = NA,
               runtime = NA,
               train_size = 0.8,
               downsample = 1.0,
               nsim = 200,
               xp = list(list(method = "radial")))
sim = bind_rows(
  sim %>% filter(dataset != "GSE_133822"),
  sim %>% filter(dataset == "GSE_133822") %>% crossing(classes = list(c("Sepsis", "Healthy"), c("Sepsis", "Crit-Ill"), c("Healthy", "Crit-Ill")))
)
sim = sim %>% arrange(dataset)

for (i in seq_len(nrow(sim))) {
  params = sim[i,]
  nsim = params$nsim
  print(i)
  print(params)

  dataset = params$dataset
  path_to_data_dir = here("processed_data", dataset)
  rds_file = list.files(path_to_data_dir, pattern = "\\.RDS$")
  print(rds_file)
  data = readRDS(here("processed_data", dataset, rds_file))
  
  if (dataset == "E_MTAB_1532") {
    predictors = data$Y
    classes = data$phenotype$Characteristics.disease.state.
  } else if (dataset == "GSE_101794") {
    predictors = data$Y
    classes = data$phenotype$`diagnosis:ch1`
  } else if (dataset == "GSE_112611") {
    BL.ix = which(data$phenotype$`baseline vs follow-up:ch1` == "BL")
    predictors = data$Y[BL.ix,]
    classes = factor(data$phenotype$`diagnosis:ch1`)[BL.ix]
  } else if (dataset == "GSE_112987") {
    predictors = data$beta
    classes = factor(data$phenotype$`disease state:ch1`)
  } else if (dataset == "GSE_85566") {
    predictors = data$beta
    classes = factor(data$phenotype$`disease status:ch1`)
  } else if (dataset == "GSE_66351") {
    predictors = data$beta
    df = data$phenotype
    classes = factor(df$`diagnosis:ch1`)
    # Define id and fctr
    id = df$`donor_id:ch1`
    ct = factor(df$`cell type:ch1`)
    br = factor(df$`brain_region:ch1`)
    fctr = sapply(1:length(ct), function(i) ifelse(ct[i] == "bulk", paste(ct[i], br[i]), paste(ct[i]) ) )
    fctr = factor(fctr)
    # Add to params for experiment_lso
    params$fctr = list(fctr)
    params$id = list(id)
  } else if (dataset == "GSE_133822") {
    predictors = data$Y
    df = data$phenotype
    # Drop "Cancer" samples (very few per cell once stratified by cell type)
    cancer_ix = which(df$patient_group=="Cancer")
    predictors = predictors[-cancer_ix,]
    df = df[-cancer_ix,]
    classes = factor(df$patient_group)
    # Subset to two classes of interest since multi-class
    ix = which(classes %in% unlist(params$classes))
    predictors = predictors[ix,]
    classes = factor(classes[ix])
    # Define id and fctr
    id = df$subject_id[ix]
    fctr = factor(df$cell_type[ix])
    # Add to params for experiment_lso
    params$fctr = list(fctr)
    params$id = list(id)
  }
  
  # Run experiment
  res = if (dataset == "GSE_66351" | dataset == "GSE_133822") {
    chunked_parLapply(nsim, experiment_lso_rt, params, predictors, classes, cl)
  } else {
    chunked_parLapply(nsim, experiment_rt, params, predictors, classes, cl) 
  }
  if (dataset == "GSE_66351" | dataset == "GSE_133822") {
    sim$subgroup_acc[[i]] = list(sapply(res, function(x) x$subgroup_acc))
  }
  sim$acc[[i]] = if (params$clf == "crc") {
    list(sapply(res, function(x) c("S" = x$accS, "L" = x$accL, "C" = x$accC)))
  } else {
    list(sapply(res, function(x) x$acc)) 
  }
  sim$runtime[[i]] = list(sapply(res, function(x) x$runtime))
  sim = as_tibble(sim)
  saveRDS(sim, here("additional_methods", "SVM_radial_66351.RDS"))
}
stopCluster(cl)