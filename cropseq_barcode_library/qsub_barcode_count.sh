#! /bin/bash
#Count barcodes from NGS reads
#Mohamad Najia

# Declare variables
project_dir=/broad/hptmp/mnajia/cropseq_barcode_library
output_dir=${project_dir}/barcode_counts
barcode_library=${project_dir}/barcode_library.tsv
samplesheet=${project_dir}/samples.txt

cd $project_dir

# Submit UGER jobs for each sample
while IFS=$'\t' read  SAMPLE  FQ1; do
  echo "Processing sample: " $SAMPLE
  qsub -N ${SAMPLE} -o $output_dir/${SAMPLE}_out.log -e $output_dir/${SAMPLE}_err.log $project_dir/barcode_count.sh $SAMPLE $FQ1 $barcode_library $output_dir
done < ${samplesheet}

