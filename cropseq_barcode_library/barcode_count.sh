#!/bin/sh
#$ -l h_vmem=15g
#$ -cwd
#$ -q broad
#$ -l h_rt=12:00:00


source /broad/software/scripts/useuse 
reuse .r-3.6.0-bioconductor
export R_LIBS_USER="/broad/blainey_lab/Mo/diffTF_env/R/3.6/libs/"

# Parse command line inputs
SAMPLE=$1
FQ1=$2
barcode_library=$3
output_dir=$4

Rscript barcode_count_v4.1.R $SAMPLE $FQ1 $barcode_library $output_dir

