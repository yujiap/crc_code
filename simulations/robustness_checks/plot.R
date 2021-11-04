library(ggplot2)
library(tibble)
library(dplyr)
library(knitr)
library(latex2exp)
sem <- function(x) return(sd(x)/sqrt(length(x)))

## Preps tibble for plotting
prep_tibble <- function(fpath) {
  sim <- readRDS(fpath)
  sim <- sim %>% select(-matches('t\\.|nsel\\.|sel'))
  
  if (grepl("other", fpath)) {
    sim$SE.acc.simple <- sapply(sim$acc.simple, function(x) x %>% unlist %>% sem)
    sim$SE.acc.uncor  <- sapply(sim$acc.uncor, function(x) x %>% unlist %>% sem)
    sim$SE.acc.corr   <- sapply(sim$acc.corr, function(x) x %>% unlist %>% sem)
    sim$acc.simple <- sapply(sim$acc.simple, function(x) x %>% unlist %>% mean)
    sim$acc.uncor <- sapply(sim$acc.uncor, function(x) x %>% unlist %>% mean)
    sim$acc.corr <- sapply(sim$acc.corr, function(x) x %>% unlist %>% mean)
  } else {
    sim$SE.acc.simple.S <- sapply(sim$acc.simple.S, function(x) x %>% unlist %>% sem)
    sim$SE.acc.simple.L <- sapply(sim$acc.simple.L, function(x) x %>% unlist %>% sem)
    sim$SE.acc.simple.C <- sapply(sim$acc.simple.C, function(x) x %>% unlist %>% sem)
    sim$SE.acc.uncor.S  <- sapply(sim$acc.uncor.S, function(x) x %>% unlist %>% sem)
    sim$SE.acc.uncor.L  <- sapply(sim$acc.uncor.L, function(x) x %>% unlist %>% sem)
    sim$SE.acc.uncor.C  <- sapply(sim$acc.uncor.C, function(x) x %>% unlist %>% sem)
    sim$SE.acc.corr.S  <- sapply(sim$acc.corr.S, function(x) x %>% unlist %>% sem)
    sim$SE.acc.corr.L  <- sapply(sim$acc.corr.L, function(x) x %>% unlist %>% sem)
    sim$SE.acc.corr.C  <- sapply(sim$acc.corr.C, function(x) x %>% unlist %>% sem)
    
    sim$acc.simple.S <- sapply(sim$acc.simple.S, function(x) x %>% unlist %>% mean)
    sim$acc.simple.L <- sapply(sim$acc.simple.L, function(x) x %>% unlist %>% mean)
    sim$acc.simple.C <- sapply(sim$acc.simple.C, function(x) x %>% unlist %>% mean)
    sim$acc.uncor.S  <- sapply(sim$acc.uncor.S, function(x) x %>% unlist %>% mean)
    sim$acc.uncor.L  <- sapply(sim$acc.uncor.L, function(x) x %>% unlist %>% mean)
    sim$acc.uncor.C  <- sapply(sim$acc.uncor.C, function(x) x %>% unlist %>% mean)
    sim$acc.corr.S  <- sapply(sim$acc.corr.S, function(x) x %>% unlist %>% mean)
    sim$acc.corr.L  <- sapply(sim$acc.corr.L, function(x) x %>% unlist %>% mean)
    sim$acc.corr.C  <- sapply(sim$acc.corr.C, function(x) x %>% unlist %>% mean) 
  }
  
  ## Prep tibble for plotting
  if (grepl("other", fpath)) {
    sim.acc <- sim %>% 
      select(-starts_with('SE.')) %>% 
      tidyr::pivot_longer(cols = matches('acc.'), 
                          names_to = 'model', 
                          values_to = 'accuracy') %>%
      mutate(model = case_when(
        grepl('simple', model) ~ 'simple',
        grepl('uncor', model) ~ 'uncorrelated',
        grepl('corr', model) ~ 'correlated'
      ))
    
    sim.SE <- sim %>% 
      select(-starts_with('acc.')) %>% 
      tidyr::pivot_longer(cols = matches('SE.'), 
                          names_to = 'model', 
                          values_to = 'SE') %>%
      mutate(model = case_when(
        grepl('simple', model) ~ 'simple',
        grepl('uncor', model) ~ 'uncorrelated',
        grepl('corr', model) ~ 'correlated'
      ))
  } else {
    sim.acc <- sim %>% 
      select(-starts_with('SE.')) %>% 
      tidyr::pivot_longer(cols = matches('acc.'), 
                          names_to = 'model', 
                          values_to = 'accuracy') %>%
      mutate(clf = case_when(
        grepl('.S', model) ~ 'S',
        grepl('.L', model) ~ 'L',
        grepl('.C', model) ~ 'C'
      )) %>%
      mutate(model = case_when(
        grepl('simple', model) ~ 'simple',
        grepl('uncor', model) ~ 'uncorrelated',
        grepl('corr', model) ~ 'correlated'
      ))
    
    sim.SE <- sim %>% 
      select(-starts_with('acc.')) %>% 
      tidyr::pivot_longer(cols = matches('SE.'), 
                          names_to = 'model', 
                          values_to = 'SE') %>%
      mutate(clf = case_when(
        grepl('.S', model) ~ 'S',
        grepl('.L', model) ~ 'L',
        grepl('.C', model) ~ 'C'
      )) %>%
      mutate(model = case_when(
        grepl('simple', model) ~ 'simple',
        grepl('uncor', model) ~ 'uncorrelated',
        grepl('corr', model) ~ 'correlated'
      ))
  }
  sim <- dplyr::left_join(sim.acc, sim.SE, by = c("n", "p", "alpha.rs.level", "L.type", "eps.type", "nsim", "model", "clf"))
  sim <- sim %>% mutate(model = factor(model, levels = c("simple", "uncorrelated", "correlated")))
  sim$alpha.rs.level <- format(sim$alpha.rs.level, scientific=FALSE)
  sim$L.type <- factor(sim$L.type)
  
  if (grepl("other", fpath)) {
    sim$clf <- factor(sim$clf, levels = c("glmnet", "pam", "myDlda"))
  } else {
    sim$clf <- factor(sim$clf, levels = c("S", "L", "C"))
  }
  return(sim)
}

## Extra line to indicate Bayes accuracy rates in all plots
extra_line = tibble(model = factor(c("Simple", "Uncorrelated", "Correlated")), yint = c(NA, NA, 0.92))

## ======================================================================================================
## L TYPE
## ======================================================================================================

## Color palette
cc = c("pam"="#cc7cf7", "C"="#287bf7", "myDlda"="#7eb215", "glmnet"="red")

## L.type: Accuracy vs n
sim_crc <- prep_tibble("sim_Ltype.RDS")
levels(sim_crc$L.type) <- c("normal", "skewnormal", "mixture", "uniform")
plt <- sim_crc %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(n != 10 & n!= 20) %>%
  ggplot(aes(x = n, y = accuracy, color = clf)) +
  facet_grid(L.type ~ model) +
  geom_line(aes(color=clf), size = 0.5, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7, size = 1, stroke = 0) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.3, width=0.3, alpha = 0.7) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  scale_color_discrete(name = "Classifier") +
  theme(legend.position="bottom") +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crc_ltype.pdf', plt, width = 8, height = 10)

sim <- prep_tibble("sim_Ltype_other.RDS")
levels(sim$L.type) <- c("normal", "skewnormal", "mixture", "uniform")
plt <- sim_crc[,order(colnames(sim_crc))] %>%
  filter(clf == "C") %>%
  rbind(sim) %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  ggplot(aes(x = n, y = accuracy, color = clf)) +
  scale_color_manual(name = "Classifier", values = cc, labels = c("myDlda" = "DLDA",
                                                                  "glmnet" = "glmnet",
                                                                  "C" = "CRC",
                                                                  "pam" = "PAM")) +
  facet_grid(L.type ~ model) +
  geom_line(aes(color=clf), size = 0.5, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7, size = 1, stroke = 0) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.3, width=0.3, alpha = 0.7) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crc_ltype_all.pdf', plt, width = 8, height = 10)

## ======================================================================================================
## ALPHA
## ======================================================================================================

## alpha: Accuracy vs n
sim_crc <- prep_tibble("sim_alpha.RDS")
# Keep only alpha = 1 result (for clarity)
sim_crc <- sim_crc[!(sim_crc$model == "simple" & as.numeric(sim_crc$alpha.rs.level) < 1),]
plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "S") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = "#F8766D", size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = "#F8766D", size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = "#F8766D", size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcS_alpha.pdf', plt, width = 10, height = 4)

plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "L") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = "#00BA38", size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = "#00BA38", size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = "#00BA38", size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcL_alpha.pdf', plt, width = 10, height = 4)

plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "C") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = "#619CFF", size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = "#619CFF", size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = "#619CFF", size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcC_alpha.pdf', plt, width = 10, height = 4)

sim <- prep_tibble("sim_alpha_other.RDS")
# Keep only alpha = 1 result (for clarity)
sim <- sim[!(sim$model == "simple" & as.numeric(sim$alpha.rs.level) < 1),]
plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "glmnet") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = cc["glmnet"], size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = cc["glmnet"], size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = cc["glmnet"], size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('glmnet_alpha.pdf', plt, width = 10, height = 4)

plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "pam") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = cc["pam"], size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = cc["pam"], size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = cc["pam"], size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('pam_alpha.pdf', plt, width = 10, height = 4)

plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "myDlda") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(alpha = alpha.rs.level), color = cc["myDlda"], size = 1, stroke = 0) +
  geom_line(aes(alpha = alpha.rs.level), color = cc["myDlda"], size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, alpha = alpha.rs.level), color = cc["myDlda"], size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(alpha = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('dlda_alpha.pdf', plt, width = 10, height = 4)

## ======================================================================================================
## ALPHA - NEW COLORS
## ======================================================================================================

# Line colors
options(scipen=999)
cc = scales::seq_gradient_pal("#EC7605", "#3310F7", "Lab")(seq(0, 1, length.out=5))
names(cc) = c("0.0001", "0.0010", "0.0100", "0.1000", "1.0000")

## alpha: Accuracy vs n
sim_crc <- prep_tibble("sim_alpha.RDS")
# Keep only alpha = 1 result (for clarity)
sim_crc <- sim_crc[!(sim_crc$model == "simple" & as.numeric(sim_crc$alpha.rs.level) < 1),]
plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "S") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) +
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcS_alpha.pdf', plt, width = 10, height = 4)

plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "L") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcL_alpha.pdf', plt, width = 10, height = 4)

plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "C") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crcC_alpha.pdf', plt, width = 10, height = 4)

sim <- prep_tibble("sim_alpha_other.RDS")
# Keep only alpha = 1 result (for clarity)
sim <- sim[!(sim$model == "simple" & as.numeric(sim$alpha.rs.level) < 1),]
plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "glmnet") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('glmnet_alpha.pdf', plt, width = 10, height = 4)

plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "pam") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('pam_alpha.pdf', plt, width = 10, height = 4)

plt <- sim %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  filter(clf == "myDlda") %>%
  ggplot(aes(x = n, y = accuracy, group = alpha.rs.level)) +
  facet_grid(.~ model) +
  geom_point(aes(color = alpha.rs.level), size = 1, stroke = 0) +
  scale_colour_manual(values = cc) +
  geom_line(aes(color = alpha.rs.level), size = 0.5) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = alpha.rs.level), size = 0.3, width=0.3) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  labs(color = TeX("$\\alpha$ sparsity level")) +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('dlda_alpha.pdf', plt, width = 10, height = 4)

## ======================================================================================================
## EPS TYPE
## ======================================================================================================

## Color palette
cc = c("pam"="#cc7cf7", "C"="#287bf7", "myDlda"="#7eb215", "glmnet"="red")

## eps.type: Accuracy vs n
sim_crc <- prep_tibble("sim_epstype.RDS")
sim_crc$eps.type <- recode(sim_crc$eps.type, dexp = "Double exponential", 
                           mixture = "Mixture", 
                           t = "t (df =3)",
                           normal = "Normal",
                           unif = "Uniform")
plt <- sim_crc %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  ggplot(aes(x = n, y = accuracy, color = clf)) +
  facet_grid(eps.type ~ model) +
  geom_line(aes(color=clf), size = 0.5, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7, size = 1, stroke = 0) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.3, width=0.3, alpha = 0.7) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  scale_color_discrete(name = "Classifier") +
  theme(legend.position="bottom") +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crc_epstype.pdf', plt, width = 8, height = 10)

sim <- prep_tibble("sim_epstype_other.RDS")
sim$eps.type <- recode(sim$eps.type, dexp = "Double exponential", 
                       mixture = "Mixture", 
                       t = "t (df =3)",
                       normal = "Normal",
                       unif = "Uniform")
plt <- sim_crc[,order(colnames(sim_crc))] %>%
  filter(clf == "C") %>%
  rbind(sim) %>%
  filter(n != 10 & n!= 20) %>%
  mutate(L.type = recode(L.type, uniform = "Uniform", normal = "Normal", mixture = "Mixture", skewnormal = "Skew normal")) %>%
  mutate(model = recode_factor(model, simple = "Simple", uncorrelated = "Uncorrelated", correlated = "Correlated")) %>%
  ggplot(aes(x = n, y = accuracy, color = clf)) +
  scale_color_manual(name = "Classifier", values = cc, labels = c("myDlda" = "DLDA",
                                                                  "glmnet" = "glmnet",
                                                                  "C" = "CRC",
                                                                  "pam" = "PAM")) +
  facet_grid(eps.type ~ model) +
  geom_line(aes(color=clf), size = 0.5, alpha = 0.7) +
  geom_point(aes(color = clf), alpha = 0.7, size = 1, stroke = 0) +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE, color = clf), size = 0.3, width=0.3, alpha = 0.7) + # size = 0.7
  geom_hline(yintercept = 0.84, color = "darkgrey", linetype="dashed", size = 0.5) +
  geom_hline(data = extra_line, aes(yintercept = yint), color = "darkgrey", linetype="dashed", size = 0.5) +
  theme_bw() +
  ylab("Accuracy") +
  theme(legend.position="bottom") +
  scale_x_sqrt() +
  ylim(0.48, 1)
ggsave('crc_epstype_all.pdf', plt, width = 8, height = 10)