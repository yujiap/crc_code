library(GEOquery)

# Get other variables
gsemat = getGEO(filename = '../../raw_data/GSE_133822/GSE133822_series_matrix.txt.gz', getGPL = FALSE)
df = pData(phenoData(gsemat))
names(df)[c(5, 8, 44, 45, 46, 47, 48)] = c('last_update_date', 'source_name', 'cell_type', 'patient_group', 'sample_id', 'sepsis_stage', 'subject_id')

df$patient_group = factor(df$patient_group)
df$cell_type = factor(df$cell_type)
df$sepsis_stage[df$sepsis_stage=="na" | df$sepsis_stage=="na_sepsis"] = NA
df$sepsis_stage = factor(df$sepsis_stage)

counts = read.table('../../raw_data/GSE_133822/GSE133822_sepsis_count.txt.gz', header = TRUE)

# Sanity check
min_counts = apply(counts[,-1], 1, min)
min(min_counts); max(min_counts)

max_counts = apply(counts[,-1], 1, max)
min(max_counts); max(max_counts)

# There are duplicated genes
length(counts[,1])
length(unique(counts[,1]))
length(counts[,1]) - length(unique(counts[,1]))

# Remove duplicated genes
ix_dup = which(table(counts[,1])>1)
length(ix_dup)
names_dup = names(table(counts[,1]))[ix_dup]
to_drop = c()
for (l in names_dup) {
    dups = which(counts[,1]==l)
    to_drop = c(to_drop, sample(dups, length(dups)-1))
}
counts_nodups = counts[-to_drop,]

# Set rownames, now that we've removed duplicated genes
rownames(counts_nodups) = counts_nodups[,1]
counts_nodups[,1] = NULL
counts_nodups[1:10, 1:3]

# Log2 transform the counts
transformData = t(log(counts_nodups + 1, base = 2))

# Account for library size
transformData_rc = as.matrix(transformData) - rowMeans(transformData)

# Save processed data
saveRDS(list("Y"=transformData_rc, "phenotype"=df), "GSE133822_processed.RDS", compress = TRUE)
