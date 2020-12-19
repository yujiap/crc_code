library(oligo)

rawset = readRDS('../../raw_data/E_MTAB_1532/E_MTAB_1532_raw.RDS')
df = pData(phenoData(rawset))
for (x in colnames(df)) {
    df[,x] = factor(df[,x])
}

# Background correction
exprMat = rma(rawset, background = TRUE, normalize = FALSE)
exprMat = t(exprs(exprMat))

# Save processed data
saveRDS(list("Y"=exprMat, "phenotype"=df), "EMTAB1532_processed.RDS", compress = TRUE)
