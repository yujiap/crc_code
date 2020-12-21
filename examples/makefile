.PHONY: GSE_85566 GSE_66351 GSE_133822 GSE_112611 GSE_112987 GSE_101794 E_MTAB_1532 examples_tables_figs

GSE_85566:
	echo "Cleaning up..."
	find ./processed_data/GSE_85566 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_85566 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_85566/ && Rscript preprocess.R
	cd ./GSE_85566/ && Rscript sim.R

GSE_66351:
	echo "Cleaning up..."
	find ./processed_data/GSE_66351 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_66351 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_66351/ && Rscript preprocess.R
	cd ./GSE_66351/ && Rscript sim_all.R

GSE_133822:
	echo "Cleaning up..."
	find ./processed_data/GSE_133822 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_133822 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_133822/ && Rscript preprocess.R
	cd ./GSE_133822/ && Rscript sim_all.R

GSE_112611:
	echo "Cleaning up..."
	find ./processed_data/GSE_112611 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_112611 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_112611/ && Rscript preprocess.R
	cd ./GSE_112611/ && Rscript sim.R

GSE_112987:
	echo "Cleaning up..."
	find ./processed_data/GSE_112987 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_112987 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_112987/ && Rscript preprocess.R
	cd ./GSE_112987/ && Rscript sim.R

GSE_101794:
	echo "Cleaning up..."
	find ./processed_data/GSE_101794 ! -name '*.R' -type f -exec rm -f {} +
	find ./GSE_101794 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/GSE_101974/ && Rscript preprocess.R
	cd ./GSE_101974/ && Rscript sim.R

E_MTAB_1532:
	echo "Cleaning up..."
	find ./processed_data/E_MTAB_1532 ! -name '*.R' -type f -exec rm -f {} +
	find ./E_MTAB_1532 ! -name '*.R' -type f -exec rm -f {} +
	# cd ./processed_data/E_MTAB_1532/ && Rscript preprocess.R
	cd ./E_MTAB_1532/ && Rscript sim.R

examples_tables_figs:
	make GSE_85566
	make GSE_66351
	make GSE_133822
	make GSE_112611
	make GSE_112987
	make GSE_101794
	make E_MTAB_1532
	Rscript table.R
	Rscript celltype_sims_table.R
	Rscript scree_plots.R
