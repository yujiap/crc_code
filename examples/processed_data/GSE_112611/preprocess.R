library(GEOquery)
library(magrittr)

# Get beta values (according to GEO, background corrected, not quantile normalized)
# Large file! Takes a few minutes to load!
zz = gzfile('../../raw_data/GSE_112611/GSE112611_beta_values.txt.gz')
dat = read.csv(zz, sep = '\t', row.names = 1, 
                   colClasses = c("character", rep('numeric', 804)))

# Remove pvalue columns
detection_pval_cols = grep("Detection.PVal*", colnames(dat))
featMat = t(dat[,-detection_pval_cols])
dim(featMat)

# Some features/observations are NA
ft_has_NA = apply(featMat, 2, function(x) sum(is.na(x))!=0)
sum(ft_has_NA)
head(ft_has_NA)
sum(is.na(featMat[2,]))
sum(is.na(featMat[3,]))
sum(is.na(featMat[,2]))
sum(is.na(featMat[,3]))

# Remove features with NA
featMat = featMat[,!ft_has_NA] # 402 504790

# Get other variables
gsemat = getGEO(filename = '../../raw_data/GSE_112611/GSE112611_series_matrix.txt.gz', getGPL = FALSE)
df = pData(phenoData(gsemat))

# Clean up: "NA" -> <NA>
df$`disease stage at diagnosis:ch1` = sapply(df$`disease stage at diagnosis:ch1`, 
                                             function(x) ifelse(x == 'NA', NA, x))

df$`disease stage at at follow up:ch1` = sapply(df$`disease stage at at follow up:ch1`, 
                                             function(x) ifelse(x == 'NA', NA, x))

df$`age:ch1` = as.numeric(df$`age:ch1`)
df$`baseline vs follow-up:ch1` = factor(df$`baseline vs follow-up:ch1`)
df$`diagnosis:ch1` = factor(df$`diagnosis:ch1`)
df$`disease stage at at follow up:ch1` = factor(df$`disease stage at at follow up:ch1`)
df$`disease stage at diagnosis:ch1` = factor(df$`disease stage at diagnosis:ch1`)
df$`gender:ch1` = factor(df$`gender:ch1`)

# Save processed data
saveRDS(list("Y"=featMat, "phenotype"=df), "GSE112611_processed.RDS", compress = TRUE)
