.PHONY: preprocess examples_tables_figs

preprocess:
	echo "Processing GSE_85566..."
	cd ./processed_data/GSE_85566/ && Rscript preprocess.R
	echo "Processing GSE_66351..."
	cd ./processed_data/GSE_66351/ && Rscript preprocess.R
	echo "Processing GSE_133822..."
	cd ./processed_data/GSE_133822/ && Rscript preprocess.R
	echo "Processing GSE_112611..."
	cd ./processed_data/GSE_112611/ && Rscript preprocess.R
	echo "Processing GSE_112987..."
	cd ./processed_data/GSE_112987/ && Rscript preprocess.R
	echo "Processing GSE_101794..."
	cd ./processed_data/GSE_101794/ && Rscript preprocess.R
	echo "Processing E_MTAB_1532..."
	cd ./processed_data/E_MTAB_1532/ && Rscript preprocess.R


examples_tables_figs:
	# make preprocess
	Rscript sim_examples.R
	Rscript table.R
	Rscript celltype_sims_table.R
	Rscript scree_plots.R