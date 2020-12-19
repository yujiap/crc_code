library(GEOquery)
library(minfi)
library(lumi)

untar('../../raw_data/GSE_85566/GSE85566_RAW.tar', exdir = './tmp_data/GSE85566_RAW/')
sapply(list.files('./tmp_data/GSE85566_RAW/', pattern = "\\.idat.gz$"), 
       function(x) gunzip(paste0('./tmp_data/GSE85566_RAW/', x)))

rgset = read.metharray.exp("./tmp_data/GSE85566_RAW/")
mset = preprocessIllumina(rgset, bg.correct=TRUE, normalize="no")
beta = getBeta(mset, type="Illumina")

gsemat = getGEO(filename = '../../raw_data/GSE_85566/GSE85566_series_matrix.txt.gz', getGPL = FALSE)
phenotype = names(pData(phenoData(gsemat)))[39:45]
df = pData(phenoData(gsemat))[, phenotype]
df$`age:ch1` = as.numeric(df$`age:ch1`)
df$`cell type:ch1` = factor(df$`cell type:ch1`)
df$`ChIP:ch1` = factor(df$`ChIP:ch1`)
df$`disease status:ch1` = factor(df$`disease status:ch1`)
df$`ethnicity:ch1` = factor(df$`ethnicity:ch1`)
df$`gender:ch1` = factor(df$`gender:ch1`)
df$`well:ch1` = factor(df$`well:ch1`)

saveRDS(list("beta"=t(beta), "phenotype"=df), "GSE85566_processed.RDS", compress=TRUE)

# Remove temp directory
unlink("./tmp_data/", recursive=TRUE)