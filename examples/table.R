library(magrittr)
library(knitr)
library(plyr)
library(dplyr)
library(stringr)

# Load data
loadRData = function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# For a given binary classification task, aggregate results across clfs
aggregateTask = function(dataset, task, clfs = c("crc", "glmnet", "pam", "myDlda")) {
  fpaths = c(list.files(path = dataset, pattern = task, full.names = TRUE))
  table.Nsim = lapply(fpaths, function(fpath) {
    results = loadRData(fpath)
    results$out
  }) %>% do.call(cbind, .)
  table.Nsim = table.Nsim[,c('glmnet', 'pam', 'myDlda', 'S', 'L', 'C')]
  table.acc = colMeans(table.Nsim)
  table.SE = apply(table.Nsim, 2, sd) / sqrt(nrow(table.Nsim))

  ret = list()
  ret$table.acc = table.acc
  ret$table.SE = table.SE
  return(ret)
}

# ---------------------------------------------------------------
# Subtables - acc
# ---------------------------------------------------------------

name_map = NULL

dataset = "GSE_133822/"
tasks = c("\\_Healthy_vs_Sepsis.out" = "Sepsis (Healthy vs. Sepsis)", 
		"\\_Crit-Ill_vs_Sepsis.out" = "Sepsis (Crit. Ill. vs. Sepsis)",
		"\\_Crit-Ill_vs_Healthy.out" = "Sepsis (Crit. Ill. vs. Healthy)")
name_map = c(name_map, tasks)
table.acc.133822 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.133822 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "E_MTAB_1532"
tasks = c(
  "\\_colorectal cancer_vs_normal.out$" = "Colorectal Cancer"
)
name_map = c(name_map, tasks)
table.acc.1352 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.1352 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "GSE_112611"
tasks = c(
  "\\_Crohn's disease_vs_non-IBD control.out$" = "Crohn's (Methylation)"
)
name_map = c(name_map, tasks)
table.acc.112611 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.112611 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "GSE_112987"
tasks = c(
  "\\_control_vs_FASD.out$" = "FASD"
)
name_map = c(name_map, tasks)
table.acc.112987 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.112987 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "GSE_101794"
tasks = c(
  "\\_CD_vs_Non-IBD.out$" = "Crohn's (Expression)"
)
name_map = c(name_map, tasks)
table.acc.101794 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.101794 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "GSE_85566"
tasks = c(
  "\\_Asthma_vs_Control.out$" = "Asthma"
)
name_map = c(name_map, tasks)
table.acc.85566 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.85566 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))


dataset = "GSE_66351/"
tasks = c(
  "\\_AD_vs_CTRL_.out$" = "Alzheimer's"
)
name_map = c(name_map, tasks)
table.acc.66351 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.acc))
table.SE.66351 = t(sapply(names(tasks), function(task) aggregateTask(dataset, task)$table.SE))

# ---------------------------------------------------------------
# Aggregate subtables
# ---------------------------------------------------------------

table.accuracies = rbind(
	  table.acc.66351,
	  table.acc.85566,
	  table.acc.1352,
	  table.acc.101794,
	  table.acc.112611,
	  table.acc.112987,
	  table.acc.133822
)

table.SE = rbind(
	  table.SE.66351,
  	  table.SE.85566,
  	  table.SE.1352,
  	  table.SE.101794,
  	  table.SE.112611,
  	  table.SE.112987,
  	  table.SE.133822
)

rownames(table.accuracies) = revalue(rownames(table.accuracies), name_map)
rownames(table.SE) = revalue(rownames(table.SE), name_map)

colnames(table.accuracies) = recode(colnames(table.accuracies), !!!c("S"="CRC-S", "L"="CRC-L", "C"="CRC", "pam"="PAM", "glmnet"="glmnet", "myDlda" = "DLDA"))
colnames(table.SE) = recode(colnames(table.SE), !!!c("S"="CRC-S", "L"="CRC-L", "C"="CRC", "pam"="PAM", "glmnet"="glmnet", "myDlda" = "DLDA"))

# Reorder cols
table.accuracies = table.accuracies[,c("glmnet", "PAM", "DLDA", "CRC", "CRC-S", "CRC-L")]
table.SE = table.SE[,c("glmnet", "PAM", "DLDA", "CRC", "CRC-S", "CRC-L")]

# To latex
kable(table.accuracies, format = "latex", digits = 2, booktabs = T) %>% cat(., file = "table_acc.txt")
kable(table.SE, format = "latex", digits = 3, booktabs = T) %>% cat(., file = "table_SE.txt")
