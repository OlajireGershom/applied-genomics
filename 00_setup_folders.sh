#!/bin/bash

mkdir environment && nano install_tools.sh

mkdir data 

cd data/

mkdir raw trimmed metadata

cd metadata/ 

nano sample_info.tsv SRR_accessions.txt

cd ../../

mkdir reference

cd reference/ 

mkdir genome annotation index

cd ../

mkdir qc

cd qc/

mkdir raw_fastqc trimmed_fastqc

cd ../

mkdir alignment

cd alignment/

mkdir sam bam logs

cd ../

mkdir counts

cd counts/

nano featurecounts_output.txt count_matrix.tsv

cd ../

mkdir scripts

cd scripts/

nano 01_download_data.sh 02_quality_control.sh 03_build_index.sh 04_alignment.sh 05_feature_counts.sh

#Add .gitkeep to empty folders so that git tracks them

find . -type d -empty -exec touch {}/.gitkeep \;

