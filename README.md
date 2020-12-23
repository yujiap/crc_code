# CRC

This repository contains the code, software, and data needed to reproduce the results in the following paper:

Pan, Y. and Gagnon-Bartsch, J. A. "Separating and Reintegrating Latent Variables to Improve Classification of Genomic Data" https://arxiv.org/abs/2012.11757

# Release notes:

## v1.0: CRC with processed data
Includes processed datasets. Each .tar.gz file should be placed in the corresponding subdirectory of the 'processed data' directory (i.e., GSE66351_processed.tar.gz in examples/processed_data/GSE_66351) then decompressed before running the code.

To reproduce both the simulations and the examples, run the following from the command line:
make all_analysis

To reproduce the simulations only, run the following from the command line:
make simulations

To reproduce the examples only, run the following from the command line:
make examples

Note that these make commands make take a long time to execute. Thus, it is recommended to run the code on a computing cluster if possible.
