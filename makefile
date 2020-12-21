.PHONY: ex crc20 crc100 crc500 pam20 pam100 pam500 dlda20 dlda100 dlda500 glmnet20 glmnet100 glmnet500 plot_sim all_analysis

ex:
	echo "Cleaning out old output..."
	find /examples/ -type f -name "table_SE.txt" -exec rm -f {} \;
	find /examples/ -type f -name "table_acc.txt" -exec rm -f {} \;
	cd /examples/ && make examples_tables

crc20:
	echo "Cleaning out old output..."
	find /simulations/sims_crc/ -type f -name "sim_p20000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_crc/ && Rscript sim_p20k.R

crc100:
	echo "Cleaning out old output..."
	find /simulations/sims_crc/ -type f -name "sim_p100000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_crc/ && Rscript sim_p100k.R

crc500:
	echo "Cleaning out old output..."
	find /simulations/sims_crc/ -type f -name "sim_p500000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_crc/ && Rscript sim_p500k.R
pam20:
	echo "Cleaning out old output..."
	find /simulations/sims_pam/ -type f -name "sim_p20000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_pam/ && Rscript sim_p20k.R

pam100:
	echo "Cleaning out old output..."
	find /simulations/sims_pam/ -type f -name "sim_p100000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_pam/ && Rscript sim_p100k.R

pam500:
	echo "Cleaning out old output..."
	find /simulations/sims_pam/ -type f -name "sim_p500000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_pam/ && Rscript sim_p500k.R

dlda20:
	echo "Cleaning out old output..."
	find /simulations/sims_dlda/ -type f -name "sim_p20000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_dlda/ && Rscript sim_p20k.R

dlda100:
	echo "Cleaning out old output..."
	find /simulations/sims_dlda/ -type f -name "sim_p100000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_dlda/ && Rscript sim_p100k.R

dlda500:
	echo "Cleaning out old output..."
	find /simulations/sims_dlda/ -type f -name "sim_p500000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_dlda/ && Rscript sim_p500k.R

glmnet20:
	echo "Cleaning out old output..."
	find /simulations/sims_glmnet/ -type f -name "sim_p20000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_glmnet/ && Rscript sim_p20k.R

glmnet100:
	echo "Cleaning out old output..."
	find /simulations/sims_glmnet/ -type f -name "sim_p100000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_glmnet/ && Rscript sim_p100k.R

glmnet500:
	echo "Cleaning out old output..."
	find /simulations/sims_glmnet/ -type f -name "sim_p500000*.RDS" -exec rm -f {} \;
	cd /simulations/sims_glmnet/ && Rscript sim_p500k.R

plot_sim:
	cd /simulations/ && Rscript plots.R
	find /simulations/ -type f -name "Rplots.pdf" -exec rm -f {} \;

all_analysis:
	make ex
	make glmnet20
	make crc20
	make dlda20
	make pam20
	make glmnet100
	make crc100
	make dlda100
	make pam100
	make glmnet500
	make crc500
	make dlda500
	make pam500
	make plot_sim
