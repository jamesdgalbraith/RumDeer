library(tidyverse)

# 
sample_list <- tibble(samples = list.files('/mnt/loki/johnston/rumdeer/james/02_read_alignment/outputs/trim_align_dedup/')) |>
    filter(endsWith(samples, suffix = "bam")) |>
    mutate(samples = sub(".bam", "", samples))

# 
family_data <-
    read_table("/mnt/loki/johnston/rumdeer/james/data/20240927_Deer_1to43_mCelEla1.fam",
               col_names = c("Family", "samples", "Father", "Mother", "Sex", "Phenotype")) |>
    dplyr::select("samples", "Sex")

males <-
    family_data |>
    filter(samples %in% sample_list$samples, Sex == 1)