library(magrittr)
library(knitr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(kableExtra)
sem = function(x) return(sd(x)/sqrt(length(x)))

sim = readRDS("sim_ts0.5.RDS")

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
  kable(format = "latex", digits = 2, booktabs = TRUE)

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
  kable(format = "latex", digits = 3, booktabs = TRUE)