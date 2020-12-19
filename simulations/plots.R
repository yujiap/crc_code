library(magrittr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(plyr)

options(scipen=999)

clf_df = function(clf_sims) {
  x_acc = sapply(1:length(clf_sims), function(i) { clf_sims[[i]]$accuracies %>% apply(., c(1,2), mean) %>% 
      `.`[,1] %>% c("n"=clf_sims[[i]]$n, .)})
  Nsim = dim(clf_sims[[1]]$accuracies)[3]
  x_SE = sapply(1:length(clf_sims), function(i) { clf_sims[[i]]$accuracies %>% {apply(., c(1,2), sd)/sqrt(Nsim)} %>% 
      `.`[,1] %>% c("n"=clf_sims[[i]]$n, .)})
  mx_acc = melt(as.data.frame(t(x_acc)), id.vars="n", variable.name = "model", value.name = "accuracy")
  mx_SE = melt(as.data.frame(t(x_SE)), id.vars = "n", variable.name = "model", value.name = "SE")
  return(merge(mx_acc, mx_SE))
}

crc_df = function(crc_sim, clf) {
  x_acc = sapply(1:length(crc_sim), function(i) { crc_sim[[i]]$accuracies %>% apply(., c(1,2), mean) %>% .[,clf] %>% c("n"=crc_sim[[i]]$n, .) })
  Nsim = dim(crc_sim[[1]]$accuracies)[3]
  x_SE = sapply(1:length(crc_sim), function(i) { crc_sim[[i]]$accuracies %>% {apply(., c(1,2), sd)/sqrt(Nsim)} %>% 
      `.`[,1] %>% c("n"=crc_sim[[i]]$n, .)})
  mx_acc = melt(as.data.frame(t(x_acc)), id.vars="n", variable.name = "model", value.name = "accuracy")
  mx_SE = melt(as.data.frame(t(x_SE)), id.vars = "n", variable.name = "model", value.name = "SE")
  return(merge(mx_acc, mx_SE))
}

dlda_sim_p20k <- list(
  sim_n50   <- readRDS("sims_dlda/sim_p20000_n50.RDS"), 
  sim_n100  <- readRDS("sims_dlda/sim_p20000_n100.RDS"),
  sim_n200  <- readRDS("sims_dlda/sim_p20000_n200.RDS"),
  sim_n500  <- readRDS("sims_dlda/sim_p20000_n500.RDS"),
  sim_n1000 <- readRDS("sims_dlda/sim_p20000_n1000.RDS")
)
dlda_sim_p100k <- list(
  sim_n50   <- readRDS("sims_dlda/sim_p1e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_dlda/sim_p1e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_dlda/sim_p1e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_dlda/sim_p1e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_dlda/sim_p1e+05_n1000.RDS")
)
dlda_sim_p500k <- list(
  sim_n50   <- readRDS("sims_dlda/sim_p5e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_dlda/sim_p5e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_dlda/sim_p5e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_dlda/sim_p5e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_dlda/sim_p5e+05_n1000.RDS")
)

pam_sim_p20k <- list(
  sim_n50  <- readRDS("sims_pam/sim_p20000_n50.RDS"),
  sim_n100 <- readRDS("sims_pam/sim_p20000_n100.RDS"),
  sim_n200 <- readRDS("sims_pam/sim_p20000_n200.RDS"),
  sim_n500 <- readRDS("sims_pam/sim_p20000_n500.RDS"),
  sim_n100 <- readRDS("sims_pam/sim_p20000_n1000.RDS")
)
pam_sim_p100k <- list(
  sim_n50   <- readRDS("sims_pam/sim_p1e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_pam/sim_p1e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_pam/sim_p1e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_pam/sim_p1e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_pam/sim_p1e+05_n1000.RDS")
)
pam_sim_p500k <- list(
  sim_n50   <- readRDS("sims_pam/sim_p5e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_pam/sim_p5e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_pam/sim_p5e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_pam/sim_p5e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_pam/sim_p5e+05_n1000.RDS")
)

crc_sim_p20k <- list(
  sim_n50   <- readRDS("sims_crc/sim_p20000_n50.RDS"),
  sim_n100  <- readRDS("sims_crc/sim_p20000_n100.RDS"),
  sim_n200  <- readRDS("sims_crc/sim_p20000_n200.RDS"),
  sim_n500  <- readRDS("sims_crc/sim_p20000_n500.RDS"),
  sim_n1000 <- readRDS("sims_crc/sim_p20000_n1000.RDS")  
)
crc_sim_p100k <- list(
  sim_n50   <- readRDS("sims_crc/sim_p1e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_crc/sim_p1e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_crc/sim_p1e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_crc/sim_p1e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_crc/sim_p1e+05_n1000.RDS")  
)
crc_sim_p500k <- list(
  sim_n50   <- readRDS("sims_crc/sim_p5e+05_n50.RDS"),
  sim_n100  <- readRDS("sims_crc/sim_p5e+05_n100.RDS"),
  sim_n200  <- readRDS("sims_crc/sim_p5e+05_n200.RDS"),
  sim_n500  <- readRDS("sims_crc/sim_p5e+05_n500.RDS"),
  sim_n1000 <- readRDS("sims_crc/sim_p5e+05_n1000.RDS")  
)

glm.gatherAlpha = function(p, alpha) {
  glm_sim = list()
  glm_sim$sim_n50   <- readRDS(paste0("sims_glmnet/sim_p", p, "_n50_alpha", alpha, ".RDS"))
  glm_sim$sim_n100  <- readRDS(paste0("sims_glmnet/sim_p", p, "_n100_alpha", alpha, ".RDS"))
  glm_sim$sim_n200  <- readRDS(paste0("sims_glmnet/sim_p", p, "_n200_alpha", alpha, ".RDS"))
  glm_sim$sim_n500  <- readRDS(paste0("sims_glmnet/sim_p", p, "_n500_alpha", alpha, ".RDS"))
  glm_sim$sim_n1000 <- readRDS(paste0("sims_glmnet/sim_p", p, "_n1000_alpha", alpha, ".RDS"))
  return(glm_sim)
}
glm_sim_p20k_0   <- glm.gatherAlpha("20000", "0")
glm_sim_p20k_0.2 <- glm.gatherAlpha("20000", "0.2")
glm_sim_p20k_0.4 <- glm.gatherAlpha("20000", "0.4")
glm_sim_p20k_0.5 <- glm.gatherAlpha("20000", "0.5")
glm_sim_p20k_0.6 <- glm.gatherAlpha("20000", "0.6")
glm_sim_p20k_0.8 <- glm.gatherAlpha("20000", "0.8")
glm_sim_p20k_1   <- glm.gatherAlpha("20000", "1")

glm_sim_p100k_0   <- glm.gatherAlpha("1e+05", "0")
glm_sim_p100k_0.2 <- glm.gatherAlpha("1e+05", "0.2")
glm_sim_p100k_0.4 <- glm.gatherAlpha("1e+05", "0.4")
glm_sim_p100k_0.5 <- glm.gatherAlpha("1e+05", "0.5")
glm_sim_p100k_0.6 <- glm.gatherAlpha("1e+05", "0.6")
glm_sim_p100k_0.8 <- glm.gatherAlpha("1e+05", "0.8")
glm_sim_p100k_1   <- glm.gatherAlpha("1e+05", "1")

glm_sim_p500k_0   <- glm.gatherAlpha("5e+05", "0")
glm_sim_p500k_0.2 <- glm.gatherAlpha("5e+05", "0.2")
glm_sim_p500k_0.4 <- glm.gatherAlpha("5e+05", "0.4")
glm_sim_p500k_0.5 <- glm.gatherAlpha("5e+05", "0.5")
glm_sim_p500k_0.6 <- glm.gatherAlpha("5e+05", "0.6")
glm_sim_p500k_0.8 <- glm.gatherAlpha("5e+05", "0.8")
glm_sim_p500k_1   <- glm.gatherAlpha("5e+05", "1")


p20k = rbind(
  cbind(crc_df(crc_sim_p20k, "C"), clf="CRC"),
  cbind(clf_df(pam_sim_p20k), clf="PAM"),
  cbind(clf_df(glm_sim_p20k_1), clf="glmnet (1.0)"),
  cbind(clf_df(glm_sim_p20k_0.8), clf="glmnet (0.8)"),
  cbind(clf_df(glm_sim_p20k_0.6), clf="glmnet (0.6)"),
  cbind(clf_df(glm_sim_p20k_0.5), clf="glmnet (0.5)"),
  cbind(clf_df(glm_sim_p20k_0.4), clf="glmnet (0.4)"),
  cbind(clf_df(glm_sim_p20k_0.2), clf="glmnet (0.2)"),
  cbind(clf_df(glm_sim_p20k_0), clf="glmnet (0.0)"),
  cbind(clf_df(dlda_sim_p20k), clf="DLDA"))
p20k$p = 20000

p100k = rbind(
  cbind(crc_df(crc_sim_p100k, "C"), clf="CRC"),
  cbind(clf_df(pam_sim_p100k), clf="PAM"),
  cbind(clf_df(glm_sim_p100k_1), clf="glmnet (1.0)"),
  cbind(clf_df(glm_sim_p100k_0.8), clf="glmnet (0.8)"),
  cbind(clf_df(glm_sim_p100k_0.6), clf="glmnet (0.6)"),
  cbind(clf_df(glm_sim_p100k_0.5), clf="glmnet (0.5)"),
  cbind(clf_df(glm_sim_p100k_0.4), clf="glmnet (0.4)"),
  cbind(clf_df(glm_sim_p100k_0.2), clf="glmnet (0.2)"),
  cbind(clf_df(glm_sim_p100k_0), clf="glmnet (0.0)"),
  cbind(clf_df(dlda_sim_p100k), clf="DLDA"))
p100k$p = 100000

p500k = rbind(
  cbind(crc_df(crc_sim_p500k, "C"), clf="CRC"),
  cbind(clf_df(pam_sim_p500k), clf="PAM"),
  cbind(clf_df(glm_sim_p500k_1), clf="glmnet (1.0)"),
  cbind(clf_df(glm_sim_p500k_0.8), clf="glmnet (0.8)"),
  cbind(clf_df(glm_sim_p500k_0.6), clf="glmnet (0.6)"),
  cbind(clf_df(glm_sim_p500k_0.5), clf="glmnet (0.5)"),
  cbind(clf_df(glm_sim_p500k_0.4), clf="glmnet (0.4)"),
  cbind(clf_df(glm_sim_p500k_0.2), clf="glmnet (0.2)"),
  cbind(clf_df(glm_sim_p500k_0), clf="glmnet (0.0)"),
  cbind(clf_df(dlda_sim_p500k), clf="DLDA"))
p500k$p = 500000

all_p = rbind(p20k, p100k, p500k)
all_p$model <- revalue(all_p$model, c("S"="Simple", "L1"="Uncorrelated", "L2"="Correlated"))
all_p$p <- revalue(as.character(all_p$p), c("20000"="p=20,000", "100000"="p=100,000", "500000"="p=500,000"))
all_p$p = factor(all_p$p, levels = c("p=20,000", "p=100,000", "p=500,000"))


# Line colors
cc = scales::seq_gradient_pal("yellow", "red", "Lab")(seq(0, 1, length.out=7))
names(cc) = paste0("glmnet (", c("0.0", "0.2", "0.4", "0.5", "0.6", "0.8", "1.0"), ")")
cc = append(cc, c("PAM"="#cc7cf7", "CRC"="#287bf7", "DLDA"="#7eb215"))

extra_line = data.frame(model = c("Simple", "Uncorrelated", "Correlated"),
                        yint = c(NA, NA, 0.92))


all_clf_grid = ggplot(all_p, aes(x=n, y=accuracy, group=clf)) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(p ~ model) +
  ylab("Accuracy") +
  scale_x_sqrt() +
  ylim(0.48, 1) + 
  scale_linetype_discrete(name = "Classifier") +
  geom_line(aes(color=clf), size=0.9, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.7, width=0.5, alpha = 0.7) +
  scale_color_manual(name = "Classifier", values = cc) +
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5)
ggsave("supplement_all_clf_grid.pdf", all_clf_grid, width = 10, height = 10)


# all clfs, p=100k, glmnet alpha=1 only
all_clfs_p100k = ggplot(subset(all_p, all_p$p=="p=100,000" & all_p$clf %in% c("CRC", "glmnet (1.0)", "DLDA", "PAM")), aes(x=n, y=accuracy, group=clf)) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(. ~ model) +
  ylab("Accuracy") +
  scale_x_sqrt() +
  ylim(0.48, 1) + 
  scale_linetype_discrete(name = "Classifier") +
  geom_line(aes(color=clf), size=0.9, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.7, width=0.5, alpha = 0.7) +
  scale_color_manual(name = "Classifier", values = cc) +
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5)
ggsave("all_clf_p100k.pdf", all_clfs_p100k, width = 10, height = 4)


# crc_grid
all_crc = rbind(
  cbind(crc_df(crc_sim_p20k, "S"), clf="CRC-S", p=20000),
  cbind(crc_df(crc_sim_p20k, "L"), clf="CRC-L", p=20000),
  cbind(crc_df(crc_sim_p20k, "C"), clf="CRC", p=20000),
  cbind(crc_df(crc_sim_p100k, "S"), clf="CRC-S", p=100000),
  cbind(crc_df(crc_sim_p100k, "L"), clf="CRC-L", p=100000),
  cbind(crc_df(crc_sim_p100k, "C"), clf="CRC", p=100000),
  cbind(crc_df(crc_sim_p500k, "S"), clf="CRC-S", p=500000),
  cbind(crc_df(crc_sim_p500k, "L"), clf="CRC-L", p=500000),
  cbind(crc_df(crc_sim_p500k, "C"), clf="CRC", p=500000)
)
all_crc$model = revalue(all_crc$model, c("S"="Simple", "L1"="Uncorrelated", "L2"="Correlated"))
all_crc$p <- revalue(as.character(all_crc$p), c("20000"="p=20,000", "100000"="p=100,000", "500000"="p=500,000"))
all_crc$p = factor(all_crc$p, levels = c("p=20,000", "p=100,000", "p=500,000"))

crc_cc = c("CRC"="#287bf7", "CRC-L"="#17B02B", "CRC-S"="#F45F5A")

crc_grid = 
  ggplot(all_crc, aes(x=n, y=accuracy, group=clf)) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(p ~ model) +
  ylab("Accuracy") +
  scale_x_sqrt() +
  ylim(0.48, 1) + 
  scale_linetype_discrete(name = "Classifier") +
  geom_line(aes(color=clf), size=0.9, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.7, width=0.5, alpha = 0.7) +
  scale_color_manual(name = "Classifier", values = crc_cc) +
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5)
ggsave("supplement_crc_grid.pdf", crc_grid, width = 10, height = 10)

# crc, p100k
crc_p100k_plt = ggplot(subset(all_crc, all_crc$p=="p=100,000"), aes(x=n, y=accuracy, group=clf)) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(.~model) +
  ylab("Accuracy") +
  scale_x_sqrt() +
  ylim(0.48, 1) + 
  scale_linetype_discrete(name = "Classifier") +
  geom_line(aes(color=clf), size=0.9, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.7, width=0.5, alpha = 0.7) +
  scale_color_manual(name = "Classifier", values = crc_cc) +
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5)
ggsave("crc_p100k.pdf", crc_p100k_plt, width = 10, height = 4)

