#! /bin/bash
#Build a custom kallisto index
#Mohamad Najia


output_dir=/broad/hptmp/mnajia/2026-01-14_CSR/github/kallisto_index
index_name=kallisto_Gag-EIF4E_custom_index.idx
cdna=$output_dir/hg19.mito.annot.cdna.Gag-EIF4E.custom
k=31

cd $output_dir
cat Gag-EIF4E.fa hg19.mito.annot.cdna > $cdna


qsub -N kallisto_index -o $output_dir/out.log -e $output_dir/err.log -m ea -M mnajia@broadinstitute.org /broad/blainey_lab/Mo/scripts/build_kallisto_index.sh $index_name $cdna $k
