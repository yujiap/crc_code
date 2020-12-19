library(GEOquery)

# Get other variables
gsemat = getGEO(filename = '../../raw_data/GSE_101794/GSE101794_series_matrix.txt.gz', getGPL = FALSE)
df = pData(phenoData(gsemat))
df$`age at diagnosis in years:ch1` = as.numeric(df$`age at diagnosis in years:ch1`)
df$`diagnosis:ch1` = factor(df$`diagnosis:ch1`)
df$`location:ch1` = factor(df$`location:ch1`)
df$`paris age:ch1` = factor(df$`paris age:ch1`)
df$`Sex:ch1` = factor(df$`Sex:ch1`)

# Counts matrix
untar('../../raw_data/GSE_101794/GSE101794_RAW.tar', exdir = './tmp_data/GSE101794_RAW/')

files = list.files('./tmp_data/GSE101794_RAW/', full.names = TRUE)

# Read GSM files
gsm_list = lapply(files, 
	function(gsm) {
		read.table(gsm, header = TRUE)
    }
)

# Check that all features match before merging
sapply(gsm_list, function(gsm) all.equal(gsm$Gene, gsm_list[[1]]$Gene))

# Merge GSM
counts = t(sapply(gsm_list, function(gsm) gsm$TPM))
colnames(counts) = gsm_list[[1]]$Gene

# Sanity check
min_counts = apply(counts, 1, min)
min(min_counts); max(min_counts)

max_counts = apply(counts, 1, max)
min(max_counts); max(max_counts)

# There are duplicated genes
length(colnames(counts))
length(unique(colnames(counts)))

# Remove duplicated genes
ix_dup = which(table(colnames(counts))>1)
to_drop = c()
for (l in names(ix_dup)) {
  dups = which(colnames(counts)==l)
  to_drop = c(to_drop, sample(dups, length(dups)-1))
}

counts_nodups = counts[,-to_drop]

# Log2 transform the counts
transformData = log(counts_nodups + 1, base = 2)

# Account for library size
transformData_rc = as.matrix(transformData) - rowMeans(transformData)

# Save processed data
saveRDS(list("Y"=transformData, "phenotype"=df), "GSE101794_processed.RDS", compress = TRUE)

# Remove temp directory
unlink("./tmp_data/", recursive=TRUE)
