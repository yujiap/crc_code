library(GEOquery)
library(minfi)
library(lumi)

untar('../../raw_data/GSE_66351/GSE66351_RAW.tar', exdir = "./tmp_data/GSE66351_RAW/")
sapply(list.files('./tmp_data/GSE66351_RAW/', pattern = "\\.idat.gz$"), 
       function(x) gunzip(paste0('./tmp_data/GSE66351_RAW/', x)))

rgset = read.metharray.exp("./tmp_data/GSE66351_RAW/")
mset = preprocessIllumina(rgset, bg.correct=TRUE, normalize="no")
beta = getBeta(mset, type="Illumina")

gsemat = getGEO(filename = "../../raw_data/GSE_66351/GSE66351_series_matrix.txt.gz", getGPL = FALSE)
df = pData(phenoData(gsemat))[,c(8,40:49)]
df$`diagnosis:ch1` = factor(df$`diagnosis:ch1`)
df$`Sex:ch1` = factor(df$`Sex:ch1`)
df$`sentrix_position:ch1` = factor(df$`sentrix_position:ch1`)
df$`sentrix_id:ch1` = factor(df$`sentrix_id:ch1`)
df$`cell type:ch1` = factor(df$`cell type:ch1`)
df$`braak_stage:ch1` = factor(df$`braak_stage:ch1`)
df$`brain_region:ch1` = factor(df$`brain_region:ch1`)
df$`age:ch1` = as.numeric(df$`age:ch1`) 

saveRDS(list("beta"=t(beta), "phenotype"=df), "GSE66351_processed.RDS", compress=TRUE)

# Remove temp directory
unlink("./tmp_data/", recursive=TRUE)