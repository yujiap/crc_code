library(tibble)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(reshape2)
sem = function(x) return(sd(x)/sqrt(length(x)))

## Read CRC results
sim = readRDS("../sim_examples.RDS")
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
accuracies = crc_acc
SE_accuracies = crc_SE

## Accuracy table - CRC
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
res0_crc = res0 %>% arrange(dataset.unique)

SE_accuracies = crc_SE

## SE table - CRC
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
res0_crc_SE = res0 %>% arrange(dataset.unique)

## Read results
res_RF = readRDS("RF.RDS") %>% mutate(xp = list(NULL)) %>% select(dataset, clf, acc, nsim, xp, classes)
res_fastkNN = readRDS("fastkNN.RDS") %>% mutate(xp = list(NULL)) %>% select(dataset, clf, acc, nsim, xp, classes)

# Read SVM results
resSVM = bind_rows(
  readRDS("SVM_linear_133822.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_101794.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_112611.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_112987.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_85566.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_66351.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_linear_1532.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_133822.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_101794.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_112611.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_112987.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_85566.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_66351.RDS") %>% select(dataset, clf, acc, nsim, xp, classes),
  readRDS("SVM_radial_1532.RDS") %>% select(dataset, clf, acc, nsim, xp, classes))

sim = rbind(res_RF, res_fastkNN, resSVM)

## Make clf names unique
sim$clf.unique = sapply(sim$xp, function(x) {
  a = unlist(x, recursive = FALSE)
  if (is.null(a)) {
    param_spec = ""
  } else {
    param_names = names(a)
    param_vals = as.character(a)
    param_spec = paste0("(", paste(paste0(param_names, "=", param_vals), collapse = ", "), ")")
  }
  return(param_spec)
})
sim$clf.unique = paste(sim$clf, sim$clf.unique)

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

## Estimate accuracy and SE
sim$SE.acc = sapply(sim$acc, function(x) x %>% unlist %>% sem %>% unlist)
sim$acc = sapply(sim$acc, function(x) x %>% unlist %>% mean %>% unlist)

## Accuracy table
res0 = sim %>% 
  select(dataset.unique, clf.unique, acc) %>% 
  pivot_wider(names_from = clf.unique, values_from = acc)

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
res0 = res0 %>% arrange(dataset.unique)

## Merge in CRC results
res0 = left_join(res0, res0_crc, by = "dataset.unique")

res = data.frame(t(res0)[-1,])
colnames(res) = res0$dataset.unique
res[,1] = as.numeric(res[,1])
res[,2] = as.numeric(res[,2])
res[,3] = as.numeric(res[,3])
res[,4] = as.numeric(res[,4])
res[,5] = as.numeric(res[,5])
res[,6] = as.numeric(res[,6])
res[,7] = as.numeric(res[,7])
res[,8] = as.numeric(res[,8])
res[,9] = as.numeric(res[,9])
res$Classifier = rownames(res)

## Shorten classifier names
x = res$Classifier
x = gsub("method=", "", x)
x = gsub("CRC-C", "CRC", x)
res$Classifier = x

## Format, to latex
rownames(res) = res$Classifier
res$Classifier = NULL
res[c("RF ", "fastkNN ", "SVM (linear)", "SVM (radial)", "CRC"),] %>%
  t %>%
  kbl(format = "latex", digits = 2, booktabs = TRUE, escape = FALSE) %>%
  row_spec(0, angle = 90) %>%
  cat(., file = "table_acc_supplement.txt")

######################################
######################################

## SE table
res0 = sim %>% 
  select(dataset.unique, clf.unique, SE.acc) %>% 
  pivot_wider(names_from = clf.unique, values_from = SE.acc)

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
res0 = res0 %>% arrange(dataset.unique)

## Merge in CRC results
res0 = left_join(res0, res0_crc_SE, by = "dataset.unique")

res = data.frame(t(res0)[-1,])
colnames(res) = res0$dataset.unique
res[,1] = as.numeric(res[,1])
res[,2] = as.numeric(res[,2])
res[,3] = as.numeric(res[,3])
res[,4] = as.numeric(res[,4])
res[,5] = as.numeric(res[,5])
res[,6] = as.numeric(res[,6])
res[,7] = as.numeric(res[,7])
res[,8] = as.numeric(res[,8])
res[,9] = as.numeric(res[,9])
res$Classifier = rownames(res)

## Shorten classifier names
x = res$Classifier
x = gsub("method=", "", x)
x = gsub("CRC-C", "CRC", x)
res$Classifier = x

## Format, to latex
rownames(res) = res$Classifier
res$Classifier = NULL
res[c("RF ", "fastkNN ", "SVM (linear)", "SVM (radial)", "CRC.SE"),] %>%
  t %>%
  kbl(format = "latex", digits = 3, booktabs = TRUE, escape = FALSE) %>%
  row_spec(0, angle = 90) %>%
  cat(., file = "table_SE_supplement.txt")
