library(ggplot2)
library(gridExtra)

screePlot = function(Z, title="") {
	Z = scale(Z, center = TRUE, scale = FALSE)
	Z.svd = svd(Z, nu=0, nv=0)
	eigvals = (Z.svd$d)^2
	pct.var = (eigvals / sum(eigvals))*100
	print(sum(pct.var[1:10]))
	df = data.frame("pc"=(1:length(pct.var)), "pct.var"=pct.var)
	ggplot(df[1:10,], aes(x=pc, y=pct.var)) + 
	scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
	ylim(c(0, 50)) +
	labs(title=title, x="", y="") +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_point()
}

data = readRDS("processed_data/GSE_101794/GSE101794_processed.RDS")
crohns.expr = screePlot(data$Y, "Crohn's (Expr.)")

data = readRDS("processed_data/GSE_112611/GSE112611_processed.RDS")
BL.ix = which(data$phenotype$`baseline vs follow-up:ch1` == "BL")
crohns.methyl = screePlot(data$Y[BL.ix,], "Crohn's (Methyl.)")

data = readRDS("processed_data/E_MTAB_1532/EMTAB1532_processed.RDS")
colorectal = screePlot(data$Y, "Colorectal Cancer")

data = readRDS("processed_data/GSE_112987/GSE112987_processed.RDS")
fasd = screePlot(data$beta, "FASD")

data = readRDS("processed_data/GSE_85566/GSE85566_processed.RDS")
asthma = screePlot(data$beta, "Asthma")

data = readRDS("processed_data/GSE_133822/GSE133822_processed.RDS")
sepsis = screePlot(data$Y, "Sepsis")
df = data$phenotype
ct = df$cell_type
sepsis.cd4 = screePlot(data$Y[ct=="CD4",], "Sepsis (CD4)")
sepsis.cd8 = screePlot(data$Y[ct=="CD8",], "Sepsis (CD8)")
sepsis.cd14 = screePlot(data$Y[ct=="CD14",], "Sepsis (CD14)")

data = readRDS("processed_data/GSE_66351/GSE66351_processed.RDS")
df = data$phenotype
ct = factor(df$`cell type:ch1`)
br = factor(df$`brain_region:ch1`)
fctr = sapply(1:length(ct), function(i) ifelse(ct[i] == "bulk", paste(ct[i], br[i]), paste(ct[i]) ) )
fctr = factor(fctr)

alz = screePlot(data$beta, "Alzheimer's")
alz.bt = screePlot(data$beta[fctr=="bulk Temporal cortex",], "Alzheimer's (BT)")
alz.bf = screePlot(data$beta[fctr=="bulk Frontal cortex",], "Alzheimer's (BF)")
alz.g = screePlot(data$beta[fctr=="Glia",], "Alzheimer's (G)")
alz.n = screePlot(data$beta[fctr=="Neuron",], "Alzheimer's (N)")

all.datasets = grid.arrange(asthma, crohns.expr, crohns.methyl, colorectal, fasd, sepsis, sepsis.cd4, sepsis.cd8, sepsis.cd14, alz, alz.bt, alz.bf, alz.g, alz.n, ncol=4)
ggsave("scree_plots_no_lines.pdf", all.datasets, width = 10, height = 10)