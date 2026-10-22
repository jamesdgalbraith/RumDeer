#!/usr/bin/bash

SAMPLE=$1

hostname
date
set -e

source /home/jgalbrai/miniforge3/etc/profile.d/conda.sh
conda activate sv_detection_env
export PATH="/home/jgalbrai/Programs/minibwa/":$PATH

# create temporary folder
SCRATCH=/scratch/$USER/$SAMPLE.trim_align_clean
mkdir -p $SCRATCH
cd $SCRATCH

# set variables (paths, threads etc.) The read directory was changed for the four different sequencing instances.
# This script requires that within the read directory ($READ_DIR), each sample has it's own subfolder containing all the paired end libraries
THREADS=20
REF_GENOME=/mnt/loki/johnston/rumdeer/james/data/GCF_910594005.1_mCerEla1.1_genomic.fna
READ_DIR=/mnt/loki/johnston/rumdeer/james/02_read_alignment/data/
OUT_DIR=/mnt/loki/johnston/rumdeer/james/02_read_alignment/outputs/trim_align_dedup/

# get actual sample name from read files. In the deer files the deer's ID follows the first underscore in the read files
SAMPLE_NAME=$(ls $READ_DIR/$SAMPLE/*_R1.fastq.gz | head -n 1 | sed 's/.*\///' | sed 's/_/ /g' | cut -f 2 -d' ')
echo $SAMPLE
full_start=`date +%s.%N`
start=$(date +%s)
### TRIM READS WITH FASTP
mkdir -p fastp/
for file in $(ls $READ_DIR/$SAMPLE/*_R1.fastq.gz)
do
    base=$(basename $file "_R1.fastq.gz")
    fastp --trim_poly_g -i ${READ_DIR}/${SAMPLE}/${base}_R1.fastq.gz -I ${READ_DIR}/${SAMPLE}/${base}_R2.fastq.gz \
                        -o fastp/${base}_R1.fastq.gz -O fastp/${base}_R2.fastq.gz
done
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo fastp took $seconds seconds


### ALIGN READS TO GENOME
mkdir -p minibwa/
start=$(date +%s)
for file in $(ls fastp/*_R1.fastq.gz)
do
    base=$(basename $file "_R1.fastq.gz")
    minibwa map -t ${THREADS} $REF_GENOME fastp/${base}_R?.fastq.gz | \
        samtools sort -o minibwa/${base}.bam -T /scratch/$USER/
    samtools index minibwa/${base}.bam 
done
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo minibwa took $seconds seconds

rm -r fastp/

### MERGE BAMS FROM EACH LANE
start=$(date +%s)
BAM_LIST=$(ls -d minibwa/*bam)
samtools merge -@ $THREADS - ${BAM_LIST} |
    samtools sort -@ ${THREADS} -o ${SAMPLE}.minibwa.bam -T /scratch/$USER/
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo merging took $seconds seconds

samtools index ${SAMPLE}.minibwa.bam
rm -r minibwa/

### DEDUPLICATE MERGED BAM
start=$(date +%s)
samtools collate -@ ${THREADS} -O -u ${SAMPLE}.minibwa.bam | \
    samtools fixmate -@ ${THREADS} -m -u - - | \
    samtools sort -@ ${THREADS} -u - | \
    samtools markdup -r -@ ${THREADS} - ${SAMPLE_NAME}.bam
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo deduplication took $seconds seconds

# remove temporary files
rm -r ${SAMPLE}.minibwa.bam

mv ${SAMPLE_NAME}.bam $OUT_DIR/${SAMPLE_NAME}.bam
samtools index $OUT_DIR/${SAMPLE_NAME}.bam

full_end=`date +%s.%N`
runtime=$( echo "$full_end - $full_start" | bc -l )
echo ${SAMPLE_NAME}" took " $runtime

conda deactivate

rm -r /mnt/loki/johnston/rumdeer/james/02_read_alignment/data/$SAMPLE/