library(knitr)
library(kableExtra)
library(reshape)

# -----------------------------
# Load simulation results

myDlda_sim <- list()
myDlda_sim$sim_n50 <- readRDS("sims_dlda/sim_p1e+05_n50.RDS")
myDlda_sim$sim_n100 <- readRDS("sims_dlda/sim_p1e+05_n100.RDS")
myDlda_sim$sim_n200 <- readRDS("sims_dlda/sim_p1e+05_n200.RDS")
myDlda_sim$sim_n500 <- readRDS("sims_dlda/sim_p1e+05_n500.RDS")
myDlda_sim$sim_n1000 <- readRDS("sims_dlda/sim_p1e+05_n1000.RDS")

pam_sim <- list()
pam_sim$sim_n50 <- readRDS("sims_pam/sim_p1e+05_n50.RDS")
pam_sim$sim_n100 <- readRDS("sims_pam/sim_p1e+05_n100.RDS")
pam_sim$sim_n200 <- readRDS("sims_pam/sim_p1e+05_n200.RDS")
pam_sim$sim_n500 <- readRDS("sims_pam/sim_p1e+05_n500.RDS")
pam_sim$sim_n1000 <- readRDS("sims_pam/sim_p1e+05_n1000.RDS")

crc_sim <- list()
crc_sim$sim_n50 <- readRDS("sims_crc/sim_p1e+05_n50.RDS")
crc_sim$sim_n100 <- readRDS("sims_crc/sim_p1e+05_n100.RDS")
crc_sim$sim_n200 <- readRDS("sims_crc/sim_p1e+05_n200.RDS")
crc_sim$sim_n500 <- readRDS("sims_crc/sim_p1e+05_n500.RDS")
crc_sim$sim_n1000 <- readRDS("sims_crc/sim_p1e+05_n1000.RDS")

glm.gatherAlpha = function(alpha) {
  glm_sim = list()
  glm_sim$sim_n50 <- readRDS(paste0("sims_glmnet/sim_p1e+05_n50_alpha", alpha, ".RDS"))
  glm_sim$sim_n100 <- readRDS(paste0("sims_glmnet/sim_p1e+05_n100_alpha", alpha, ".RDS"))
  glm_sim$sim_n200 <- readRDS(paste0("sims_glmnet/sim_p1e+05_n200_alpha", alpha, ".RDS"))
  glm_sim$sim_n500 <- readRDS(paste0("sims_glmnet/sim_p1e+05_n500_alpha", alpha, ".RDS"))
  glm_sim$sim_n1000 <- readRDS(paste0("sims_glmnet/sim_p1e+05_n1000_alpha", alpha, ".RDS"))
  return(glm_sim)
}

glm_sim_0.2 <- glm.gatherAlpha("0.2")
glm_sim_0.4 <- glm.gatherAlpha("0.4")
glm_sim_0.5 <- glm.gatherAlpha("0.5")
glm_sim_0.6 <- glm.gatherAlpha("0.6")
glm_sim_0.8 <- glm.gatherAlpha("0.8")
glm_sim_1 <- glm.gatherAlpha("1")

# -----------------------------
# Aggregate info about median # genes selected

Nsel.table = rbind(
  sapply(glm_sim_1, function(x) x$accuracies %>% apply(., c(1,2), median) %>% .[,"Nsel"]),
  sapply(pam_sim, function(x) x$accuracies %>% apply(., c(1,2), median) %>% .[,"Nsel"]),
  sapply(myDlda_sim, function(x) x$accuracies %>% apply(., c(1,2), median) %>% .[,"Nsel"]),
  sapply(crc_sim, function(x) x$accuracies %>% apply(., c(1,2), median) %>% .[,"Nsel"])
) %>% t
colnames(Nsel.table) = rep(c("S", "U", "C"), 4)
rownames(Nsel.table) = c("$n=50$", "$n=100$", "$n=200$", "$n=500$", "$n=1000$")


kable(Nsel.table, "latex", booktabs = TRUE, escape = FALSE) %>% 
  t %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "glmnet ($\alpha = 1$)" = 3, "PAM" = 3, "DLDA" = 3, "CRC" = 3), escape = FALSE)

