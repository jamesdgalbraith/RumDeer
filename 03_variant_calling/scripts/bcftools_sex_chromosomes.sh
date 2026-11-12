#!/usr/bin/bash

#SBATCH --job-name=bcftools_x    # Job name
#SBATCH --export=ALL
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=20            # CPU request
#SBATCH --mem=64G                     # Job memory request
#SBATCH --time=7-00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/bcftools_x.%j.out

hostname
date
set -e

## set variables (sample (chromosome to call) from array script)
SAMPLE=chr_X
THREADS=20
GENOME="/mnt/loki/johnston/rumdeer/james/data/GCF_910594005.1_mCerEla1.1_genomic.fna"
BED_DIR=/mnt/loki/johnston/rumdeer/james/03_variant_calling/data/beds/
BAM_DIR="/mnt/loki/johnston/rumdeer/james/02_read_alignment/outputs/trim_align_dedup/"
DATA_DIR="/mnt/loki/johnston/rumdeer/james/data/"
OUT_DIR="/mnt/loki/johnston/rumdeer/james/03_variant_calling/outputs/bcftools/"

echo $SAMPLE

# activate conda environment and make scrtach dir
source /home/jgalbrai/miniforge3/etc/profile.d/conda.sh
conda activate sv_detection_env

SCRATCH=/scratch/$USER/deer_bcftools_${SAMPLE}
mkdir -p $SCRATCH/bams/
cd $SCRATCH

# extract chromosome alignments for each sample
echo "Starting "$SAMPLE" samtools"
parallel -j ${THREADS} samtools view -L ${BED_DIR}/${SAMPLE}.bed \
    ${BAM_DIR}/{}.bam \
    -o bams/{}.bam --write-index ::: $(ls ${BAM_DIR}/*bam | sed 's/.*\///;s/.bam//')

# enter sampled bam dir and call variants
# ploidy and sample sex files only needed for calling sex chromosomes and mitochondria
cd bams/
ls *bam > bam_list.txt
echo "Starting "$SAMPLE" bcftools"
bcftools mpileup --threads ${THREADS} -Ou -f ${GENOME} -g 0,5,10,20,40 -a DP --max-depth 100 --min-MQ 30 \
    -b bam_list.txt -R ${BED_DIR}/${SAMPLE}.bed | \
    bcftools call --threads ${THREADS} -m -g 0,5,10,20,40 \
    --ploidy-file ${DATA_DIR}/ploidy_file.txt \
    --samples-file ${DATA_DIR}/sample_sex.tsv \
    -Oz -o ${OUT_DIR}/${SAMPLE}.mpileup.gvcf.gz --write-index

rm -r $SCRATCH

cd $HOME


