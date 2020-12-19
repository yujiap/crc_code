library(knitr)
library(kableExtra)
library(reshape)

computeSE = function(results.out) {
  apply(results.out, 2, function(x) sd(x)/sqrt(nrow(results.out)))
}

# Alzheimers

load("./GSE_66351/crc_AD_vs_CTRL_.out")
alz.crc = colMeans(results$out)[c("C.bulk Frontal cortex", "C.bulk Temporal cortex", "C.Glia", "C.Neuron")]
alz.crc = rename(alz.crc, c("C.bulk Frontal cortex"="BF", 
                  "C.bulk Temporal cortex"="BT",
                  "C.Glia"="G",
                  "C.Neuron"="N"))
load("./GSE_66351/glmnet_AD_vs_CTRL_.out")
alz.glmnet = colMeans(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

load("./GSE_66351/pam_AD_vs_CTRL_.out")
alz.pam = colMeans(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

load("./GSE_66351/myDlda_AD_vs_CTRL_.out")
alz.dlda = colMeans(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

name_map = c("bulk Frontal cortex"="BF", 
             "bulk Temporal cortex"="BT",
             "Glia"="G",
             "Neuron"="N")
alz.glmnet = rename(alz.glmnet, name_map)
alz.pam = rename(alz.pam, name_map)
alz.dlda = rename(alz.dlda, name_map)

alz.table = c(alz.glmnet, alz.pam, alz.dlda, alz.crc)
alz.names = names(alz.table)
alz.table = matrix(alz.table, nrow=1)
colnames(alz.table) = alz.names
kable(alz.table, "latex", booktabs = TRUE, escape = FALSE, digits = 2) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c("glmnet ($\alpha = 1$)" = 4, "PAM" = 4, "DLDA" = 4, "CRC" = 4), escape = FALSE)

# Alzheimer's SE

load("./GSE_66351/crc_AD_vs_CTRL_.out")
alz.crc = computeSE(results$out)[c("C.bulk Frontal cortex", "C.bulk Temporal cortex", "C.Glia", "C.Neuron")]
alz.crc = rename(alz.crc, c("C.bulk Frontal cortex"="BF", 
                            "C.bulk Temporal cortex"="BT",
                            "C.Glia"="G",
                            "C.Neuron"="N"))
load("./GSE_66351/glmnet_AD_vs_CTRL_.out")
alz.glmnet = computeSE(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

load("./GSE_66351/pam_AD_vs_CTRL_.out")
alz.pam = computeSE(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

load("./GSE_66351/myDlda_AD_vs_CTRL_.out")
alz.dlda = computeSE(results$out)[c("bulk Frontal cortex", "bulk Temporal cortex", "Glia", "Neuron")]

name_map = c("bulk Frontal cortex"="BF", 
             "bulk Temporal cortex"="BT",
             "Glia"="G",
             "Neuron"="N")
alz.glmnet = rename(alz.glmnet, name_map)
alz.pam = rename(alz.pam, name_map)
alz.dlda = rename(alz.dlda, name_map)

alz.table = c(alz.glmnet, alz.pam, alz.dlda, alz.crc)
alz.names = names(alz.table)
alz.table = matrix(alz.table, nrow=1)
colnames(alz.table) = alz.names
kable(alz.table, "latex", booktabs = TRUE, escape = FALSE, digits = 3) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c("glmnet ($\alpha = 1$)" = 4, "PAM" = 4, "DLDA" = 4, "CRC" = 4), escape = FALSE)


# Sepsis

tasks = c("Healthy_vs_Sepsis.out", "Crit-Ill_vs_Sepsis.out", "Crit-Ill_vs_Healthy.out")
sepsis.table = sapply(tasks, function(task) {
  load(paste0("GSE_133822/crc_", task))
  sepsis.crc = colMeans(results$out)[c("C.CD14", "C.CD4", "C.CD8")]
  sepsis.crc = rename(sepsis.crc, c("C.CD14"="CD14", 
                                    "C.CD4"="CD4",
                                    "C.CD8"="CD8"))
  
  load(paste0("GSE_133822/glmnet_", task))
  sepsis.glmnet = colMeans(results$out)[c("CD14", "CD4", "CD8")]
  
  load(paste0("GSE_133822/pam_", task))
  sepsis.pam = colMeans(results$out)[c("CD14", "CD4", "CD8")]
  
  load(paste0("GSE_133822/myDlda_", task))
  sepsis.dlda = colMeans(results$out)[c("CD14", "CD4", "CD8")]
  
  c(sepsis.glmnet, sepsis.pam, sepsis.dlda, sepsis.crc)
}) %>% t


kable(sepsis.table, "latex", booktabs = TRUE, escape = FALSE, digits = 2) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c(" "=1, "glmnet ($\alpha = 1$)" = 3, "PAM" = 3, "DLDA" = 3, "CRC" = 3), escape = FALSE)

# Sepsis SE

tasks = c("Healthy_vs_Sepsis.out", "Crit-Ill_vs_Sepsis.out", "Crit-Ill_vs_Healthy.out")
sapply(tasks, function(task) {

  load(paste0("GSE_133822/glmnet_", task))
  sepsis.glmnet = computeSE(results$out)[c("CD14", "CD4", "CD8")]
  
  load(paste0("GSE_133822/pam_", task))
  sepsis.pam = computeSE(results$out)[c("CD14", "CD4", "CD8")]

  c(sepsis.glmnet, sepsis.pam)
}) %>% t %>% kable(., "latex", booktabs = TRUE, escape = FALSE, digits = 3) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c(" "=1, "glmnet ($\alpha = 1$)" = 3, "PAM" = 3), escape = FALSE)

sapply(tasks, function(task) {
  
  load(paste0("GSE_133822/crc_", task))
  sepsis.crc = computeSE(results$out)[c("C.CD14", "C.CD4", "C.CD8")]
  sepsis.crc = rename(sepsis.crc, c("C.CD14"="CD14", 
                                    "C.CD4"="CD4",
                                    "C.CD8"="CD8"))
  
  load(paste0("GSE_133822/myDlda_", task))
  sepsis.dlda = computeSE(results$out)[c("CD14", "CD4", "CD8")]
  
  c(sepsis.dlda, sepsis.crc)
}) %>% t %>% kable(., "latex", booktabs = TRUE, escape = FALSE, digits = 3) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c(" "=1, "DLDA" = 3, "CRC" = 3), escape = FALSE)