library(GEOquery)
library(minfi)
library(lumi)

untar('../../raw_data/GSE_112987/GSE112987_RAW.tar', exdir = './tmp_data/GSE112987_RAW/')
sapply(list.files('./tmp_data/GSE112987_RAW/', pattern = "\\.idat.gz$"), 
       function(x) gunzip(paste0('./tmp_data/GSE112987_RAW/', x)))

rgset = read.metharray.exp("./tmp_data/GSE112987_RAW/")
mset = preprocessIllumina(rgset, bg.correct=TRUE, normalize="no")
beta = getBeta(mset, type="Illumina")

gsemat = getGEO(filename = '../../raw_data/GSE_112987/GSE112987_series_matrix.txt.gz', getGPL = FALSE)
phenotype = names(pData(phenoData(gsemat)))
df = pData(phenoData(gsemat))[, phenotype]

df$`facial:ch1` = sapply(df$`facial:ch1`, function(x) ifelse(x=='NA', NA, x))
df$`growth:ch1` = sapply(df$`growth:ch1`, function(x) ifelse(x=='NA', NA, x))
df$`cns:ch1` = sapply(df$`cns:ch1`, function(x) ifelse(x=='NA', NA, x))

df$`age:ch1` = as.numeric(df$`age:ch1`)
df$`disease state:ch1` = factor(df$`disease state:ch1`)
df$`cns:ch1` = factor(df$`cns:ch1`)
df$`facial:ch1` = factor(df$`facial:ch1`)
df$`gender:ch1` = factor(df$`gender:ch1`)
df$`growth:ch1` = factor(df$`growth:ch1`)



saveRDS(list("beta"=t(beta), "phenotype"=df), "GSE112987_processed.RDS", compress=TRUE)

# Remove temp directory
unlink("./tmp_data/", recursive=TRUE)