#! /bin/bash
#$ -cwd
#$ -q broad
#$ -l h_rt=12:00:00
#$ -l h_vmem=50g


source /broad/software/scripts/useuse
use STAR

project_dir=/broad/hptmp/mnajia/hg38_STAR_index

cd $project_dir

#Download the hg38 fasta and transcriptome GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gunzip gencode.v47.annotation.gtf.gz

fa_file=${project_dir}/GRCh38.primary_assembly.genome.fa
gtf_file=${project_dir}/gencode.v47.annotation.gtf

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ${project_dir} \
  --genomeFastaFiles ${fa_file} \
  --sjdbGTFfile ${gtf_file} \
  --sjdbOverhang 149

