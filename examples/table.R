library(magrittr)
library(knitr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(kableExtra)
sem = function(x) return(sd(x)/sqrt(length(x)))

sim = readRDS("sim_examples.RDS")

## Make dataset names unique
sim$dataset.unique = sapply(sim$classes, function(x) {
  a = unlist(x, recursive = FALSE)
  if (is.null(a)) {
    classes_spec = ""
  } else {
    classes_spec = paste0("(", a[1], " vs. ", a[2],")")
  }
  return(classes_spec)
})
sim$dataset.unique = paste(sim$dataset, sim$dataset.unique)

## Estimate accuracy and SE for CRC
sim_crc = sim %>% filter(clf == "crc")
sim_crc$S.SE = sapply(sim_crc$acc, function(x) x[[1]]["S",] %>% sd / sqrt(length(x[[1]]["S",])))
sim_crc$L.SE = sapply(sim_crc$acc, function(x) x[[1]]["L",] %>% sd / sqrt(length(x[[1]]["L",])))
sim_crc$C.SE = sapply(sim_crc$acc, function(x) x[[1]]["C",] %>% sd / sqrt(length(x[[1]]["C",])))
sim_crc$S = sapply(sim_crc$acc, function(x) x[[1]]["S",] %>% mean)
sim_crc$L = sapply(sim_crc$acc, function(x) x[[1]]["L",] %>% mean)
sim_crc$C = sapply(sim_crc$acc, function(x) x[[1]]["C",] %>% mean)

crc_acc = sim_crc %>% 
  select(dataset.unique, clf, S, L, C) %>% 
  melt(value.name = "acc", variable.name = "type") %>%
  mutate(clf = paste0("CRC-", type)) %>%
  select(-type)

crc_SE = sim_crc %>%
  select(dataset.unique, clf, S.SE, L.SE, C.SE) %>% 
  melt(value.name = "SE.acc", variable.name = "type") %>%
  mutate(clf = paste0("CRC-", type)) %>%
  select(-type)

## Estimate accuracy and SE for other classifiers
sim_other = sim %>% filter(clf != "crc")
sim_other$SE.acc = sapply(sim_other$acc, function(x) x %>% unlist %>% sem %>% unlist)
sim_other$acc = sapply(sim_other$acc, function(x) x %>% unlist %>% mean %>% unlist)

## Merge results
accuracies = bind_rows(sim_other %>% select(dataset.unique, clf, acc), crc_acc)
SE_accuracies = bind_rows(sim_other %>% select(dataset.unique, clf, SE.acc), crc_SE)

## Accuracy table
res0 = accuracies %>% pivot_wider(names_from = clf, values_from = acc)

res0$dataset.unique = recode(res0$dataset.unique,
                             "E_MTAB_1532 " = "Colorectal Cancer",
                             "GSE_101794 " = "Crohn's (Expr.)",
                             "GSE_112611 " = "Crohn's (Methyl.)",
                             "GSE_112987 " = "FASD",
                             "GSE_85566 " = "Asthma",
                             "GSE_66351 " = "Alzheimer's",
                             "GSE_133822 (Sepsis vs. Healthy)" = "Sepsis (Sepsis vs. Healthy)",
                             "GSE_133822 (Sepsis vs. Crit-Ill)" = "Sepsis (Sepsis vs. Crit-Ill)",
                             "GSE_133822 (Healthy vs. Crit-Ill)" = "Sepsis (Healthy vs. Crit-Ill)")

colnames(res0) = recode(colnames(res0),
                        "dataset.unique" = "Dataset",
                        "pam" = "PAM",
                        "CRC-C" = "CRC",
                        "myDlda" = "DLDA")

res0 = res0 %>% select(Dataset, glmnet, PAM, DLDA, `CRC`, `CRC-S`, `CRC-L`)
res0 %>% 
  arrange(Dataset) %>%
  kable(format = "latex", digits = 2, booktabs = TRUE) %>% cat(., file = "table_acc.txt")

## SE table
res0 = SE_accuracies %>% pivot_wider(names_from = clf, values_from = SE.acc)

res0$dataset.unique = recode(res0$dataset.unique,
                             "E_MTAB_1532 " = "Colorectal Cancer",
                             "GSE_101794 " = "Crohn's (Expr.)",
                             "GSE_112611 " = "Crohn's (Methyl.)",
                             "GSE_112987 " = "FASD",
                             "GSE_85566 " = "Asthma",
                             "GSE_66351 " = "Alzheimer's",
                             "GSE_133822 (Sepsis vs. Healthy)" = "Sepsis (Sepsis vs. Healthy)",
                             "GSE_133822 (Sepsis vs. Crit-Ill)" = "Sepsis (Sepsis vs. Crit-Ill)",
                             "GSE_133822 (Healthy vs. Crit-Ill)" = "Sepsis (Healthy vs. Crit-Ill)")

colnames(res0) = recode(colnames(res0),
                        "dataset.unique" = "Dataset",
                        "pam" = "PAM",
                        "CRC-C.SE" = "CRC",
                        "CRC-S.SE" = "CRC-S",
                        "CRC-L.SE" = "CRC-L",
                        "myDlda" = "DLDA")

res0 = res0 %>% select(Dataset, glmnet, PAM, DLDA, `CRC`, `CRC-S`, `CRC-L`)
res0 %>% 
  arrange(Dataset) %>%
  kable(format = "latex", digits = 3, booktabs = TRUE) %>% cat(., file = "table_SE_acc.txt")

## Cell type simulation tables

sim_crc$sg_SE = sapply(sim_crc$subgroup_acc, function(u) {
  x = u[[1]]
  if (!is.null(dim(x))) apply(x, 1, sd) / sqrt(ncol(x))
})
sim_crc$sg_mean_acc = sapply(sim_crc$subgroup_acc, function(u) {
  x = u[[1]]
  if (!is.null(dim(x))) rowMeans(x)
})
sim_other$sg_SE = sapply(sim_other$subgroup_acc, function(u) {
  x = u[[1]]
  if (!is.null(dim(x))) apply(x, 1, function(z) sd(z, na.rm = TRUE)) / sqrt(ncol(x))
})
sim_other$sg_mean_acc = lapply(sim_other$subgroup_acc, function(u) {
  x = u[[1]]
  if (!is.null(dim(x))) rowMeans(x, na.rm = TRUE)
})

## Alzheimer's accuracy
## NOTE: For PAM, one iteration returns NAs for two cell types
# x = filter(sim_other, dataset == "GSE_66351")[2,]$subgroup_acc[[1]][[1]]
# sum(is.na(x)) # 2
# x[,985] # NA for Glia and Neuron

alz = bind_rows(filter(sim_crc, dataset == "GSE_66351") %>% select(dataset.unique, clf, sg_mean_acc),
                filter(sim_other, dataset == "GSE_66351") %>% select(dataset.unique, clf, sg_mean_acc))
alz_acc = sapply(alz$sg_mean_acc, unlist)
colnames(alz_acc) = alz$clf
alz_acc = t(alz_acc)
colnames(alz_acc) = colnames(alz_acc) %>% 
  recode("bulk Frontal cortex"="BF", 
         "bulk Temporal cortex"="BT",
         "Glia"="G",
         "Neuron"="N")
alz_table = c(alz_acc["glmnet",], alz_acc["pam",], alz_acc["myDlda",], alz_acc["crc",])
alz_names = names(alz_table)
alz_table = matrix(alz_table, nrow = 1)
colnames(alz_table) = alz_names

kable(alz_table, "latex", booktabs = TRUE, escape = FALSE, digits = 2) %>%
  kable_styling() %>%
  add_header_above(c("glmnet ($\alpha = 1$)" = 4, "PAM" = 4, "DLDA" = 4, "CRC" = 4), escape = FALSE)

## Alzheimer's SE
alz = bind_rows(filter(sim_crc, dataset == "GSE_66351") %>% select(dataset.unique, clf, sg_SE),
                filter(sim_other, dataset == "GSE_66351") %>% select(dataset.unique, clf, sg_SE))
alz_SE = sapply(alz$sg_SE, unlist)
colnames(alz_SE) = alz$clf
alz_SE = t(alz_SE)
colnames(alz_SE) = colnames(alz_SE) %>% 
  recode("bulk Frontal cortex"="BF", 
         "bulk Temporal cortex"="BT",
         "Glia"="G",
         "Neuron"="N")
alz_table = c(alz_SE["glmnet",], alz_SE["pam",], alz_SE["myDlda",], alz_SE["crc",])
alz_names = names(alz_table)
alz_table = matrix(alz_table, nrow = 1)
colnames(alz_table) = alz_names
kable(alz_table, "latex", booktabs = TRUE, escape = FALSE, digits = 3) %>%
  kable_styling() %>%
  add_header_above(c("glmnet ($\alpha = 1$)" = 4, "PAM" = 4, "DLDA" = 4, "CRC" = 4), escape = FALSE)

## Sepsis accuracy
sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Sepsis vs. Healthy)") %>% select(dataset.unique, clf, sg_mean_acc),
                   filter(sim_other, dataset.unique == "GSE_133822 (Sepsis vs. Healthy)") %>% select(dataset.unique, clf, sg_mean_acc))
sepsis_acc = sapply(sepsis$sg_mean_acc, unlist)
colnames(sepsis_acc) = sepsis$clf
sepsis_acc = t(sepsis_acc)
sepsis_table = c(sepsis_acc["glmnet",], sepsis_acc["pam",], sepsis_acc["myDlda",], sepsis_acc["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table1 = sepsis_table
rownames(sepsis_table1) = "(Sepsis vs. Healthy)"

sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Sepsis vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_mean_acc),
                   filter(sim_other, dataset.unique == "GSE_133822 (Sepsis vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_mean_acc))
sepsis_acc = sapply(sepsis$sg_mean_acc, unlist)
colnames(sepsis_acc) = sepsis$clf
sepsis_acc = t(sepsis_acc)
sepsis_table = c(sepsis_acc["glmnet",], sepsis_acc["pam",], sepsis_acc["myDlda",], sepsis_acc["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table2 = sepsis_table
rownames(sepsis_table2) = "(Sepsis vs. Crit-Ill)"

sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Healthy vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_mean_acc),
                   filter(sim_other, dataset.unique == "GSE_133822 (Healthy vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_mean_acc))
sepsis_acc = sapply(sepsis$sg_mean_acc, unlist)
colnames(sepsis_acc) = sepsis$clf
sepsis_acc = t(sepsis_acc)
sepsis_table = c(sepsis_acc["glmnet",], sepsis_acc["pam",], sepsis_acc["myDlda",], sepsis_acc["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table3 = sepsis_table
rownames(sepsis_table3) = "(Healthy vs. Crit-Ill)"

sepsis_table = rbind(sepsis_table1, sepsis_table2, sepsis_table3)
kable(sepsis_table, "latex", booktabs = TRUE, escape = FALSE, digits = 2, row.names = TRUE) %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "glmnet ($\alpha = 1$)" = 3, "PAM" = 3, "DLDA" = 3, "CRC" = 3), escape = FALSE)

## Sepsis SE
sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Sepsis vs. Healthy)") %>% select(dataset.unique, clf, sg_SE),
                   filter(sim_other, dataset.unique == "GSE_133822 (Sepsis vs. Healthy)") %>% select(dataset.unique, clf, sg_SE))
sepsis_SE = sapply(sepsis$sg_SE, unlist)
colnames(sepsis_SE) = sepsis$clf
sepsis_SE = t(sepsis_SE)
sepsis_table = c(sepsis_SE["glmnet",], sepsis_SE["pam",], sepsis_SE["myDlda",], sepsis_SE["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table1 = sepsis_table
rownames(sepsis_table1) = "(Sepsis vs. Healthy)"

sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Sepsis vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_SE),
                   filter(sim_other, dataset.unique == "GSE_133822 (Sepsis vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_SE))
sepsis_SE = sapply(sepsis$sg_SE, unlist)
colnames(sepsis_SE) = sepsis$clf
sepsis_SE = t(sepsis_SE)
sepsis_table = c(sepsis_SE["glmnet",], sepsis_SE["pam",], sepsis_SE["myDlda",], sepsis_SE["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table2 = sepsis_table
rownames(sepsis_table2) = "(Sepsis vs. Crit-Ill)"

sepsis = bind_rows(filter(sim_crc, dataset.unique == "GSE_133822 (Healthy vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_SE),
                   filter(sim_other, dataset.unique == "GSE_133822 (Healthy vs. Crit-Ill)") %>% select(dataset.unique, clf, sg_SE))
sepsis_SE = sapply(sepsis$sg_SE, unlist)
colnames(sepsis_SE) = sepsis$clf
sepsis_SE = t(sepsis_SE)
sepsis_table = c(sepsis_SE["glmnet",], sepsis_SE["pam",], sepsis_SE["myDlda",], sepsis_SE["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
colnames(sepsis_table) = sepsis_names
sepsis_table3 = sepsis_table
rownames(sepsis_table3) = "(Healthy vs. Crit-Ill)"

sepsis_table = rbind(sepsis_table1, sepsis_table2, sepsis_table3)
kable(sepsis_table, "latex", booktabs = TRUE, escape = FALSE, digits = 3, row.names = TRUE) %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "glmnet ($\alpha = 1$)" = 3, "PAM" = 3, "DLDA" = 3, "CRC" = 3), escape = FALSE)








# # Sepsis SE
# sepsis = bind_rows(filter(sim_crc, dataset == "GSE_133822") %>% select(dataset.unique, clf, sg_SE),
#                    filter(sim_other, dataset == "GSE_133822") %>% select(dataset.unique, clf, sg_SE))
sepsis_SE = sapply(sepsis$sg_SE, unlist)
colnames(sepsis_SE) = sepsis$clf
sepsis_SE = t(sepsis_SE)
sepsis_table = c(sepsis_SE["glmnet",], sepsis_SE["pam",], sepsis_SE["myDlda",], sepsis_SE["crc",])
sepsis_names = names(sepsis_table)
sepsis_table = matrix(sepsis_table, nrow = 1)
# colnames(sepsis_table) = sepsis_names
# kable(sepsis_table, "latex", booktabs = TRUE, escape = FALSE, digits = 3) %>%
#   kable_styling() %>%
#   add_header_above(c("glmnet ($\alpha = 1$)" = 3, "PAM" = 3, "DLDA" = 3, "CRC" = 3), escape = FALSE)

