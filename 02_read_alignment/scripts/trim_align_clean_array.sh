#!/usr/bin/bash

#SBATCH --job-name=trim_align_clean    # Job name
#SBATCH --export=ALL
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=20            # CPU request
#SBATCH --mem=32G                     # Job memory request
#SBATCH --time=7-00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/trim_align_clean.%A.%a.out
#SBATCH --array=1-79%5

Param=`cat /mnt/loki/johnston/rumdeer/james/02_read_alignment/data/sample_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}"`
bash /mnt/loki/johnston/rumdeer/james/02_read_alignment/scripts/trim_align_clean.sh $Param