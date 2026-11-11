conda activate sv_detect_2

# get X chromosomes for males
for SAMPLE in ALD99 ART00 AST91 ATH91 BSY02; do
    echo $SAMPLE
    samtools view /mnt/loki/johnston/rumdeer/james/02_read_alignment/outputs/trim_align_dedup/${SAMPLE}.bam \
        -o /mnt/loki/johnston/rumdeer/james/01_par_check/outputs/male_X_bams/${SAMPLE}.bam --write-index -b NC_057848.1
done
