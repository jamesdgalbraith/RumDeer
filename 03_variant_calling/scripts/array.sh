#!/usr/bin/bash

#SBATCH --job-name=bcftools    # Job name
#SBATCH --export=ALL
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=20            # CPU request
#SBATCH --mem=64G                     # Job memory request
#SBATCH --time=7-00:00               # Time limit hrs:min:sec
#SBATCH --output=logs/bcftools.%A.%a.out
#SBATCH --array=1-33%3

source /home/jgalbrai/miniforge3/etc/profile.d/conda.sh
conda activate sv_detection_env

Param=`cat /mnt/loki/johnston/rumdeer/james/03_variant_calling/data/beds.txt | awk "NR==${SLURM_ARRAY_TASK_ID}"`
bash /mnt/loki/johnston/rumdeer/james/03_variant_calling/scripts/bcftools_autosomes.sh $Param