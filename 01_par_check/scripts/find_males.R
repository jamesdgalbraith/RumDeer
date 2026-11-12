library(tidyverse)

# read in pedigree, determine sex based on whether mother or father, filter out those lacking data
pedigree <- read_csv("/mnt/loki/johnston/rumdeer/james/data/Pedigree_2024-06-06_DB.csv") |>
    dplyr::mutate(PedigreeSex = case_when(Code %in% Mother ~ "F",
                                          Code %in% Father ~ "M",
                                          TRUE ~ NA)) |>
    dplyr::select(Code, PedigreeSex) |>
    dplyr::filter(!is.na(PedigreeSex)) |>
    dplyr::rename(samples = Code)
 
#  get list of samples from alignments
sample_list <- tibble(samples = list.files('/mnt/loki/johnston/rumdeer/james/02_read_alignment/outputs/trim_align_dedup/')) |>
    filter(endsWith(samples, suffix = "bam")) |>
    mutate(samples = sub(".bam", "", samples))

# read in PLINK .fam data
family_data <-
    read_table("/mnt/loki/johnston/rumdeer/james/data/20240927_Deer_1to43_mCelEla1.fam",
               col_names = c("Family", "samples", "Father", "Mother", "Sex", "Phenotype")) |>
    dplyr::select("samples", "Sex")

# get sequenced samples sex from PLINK data
samples_sex <-
    family_data |>
    filter(samples %in% sample_list$samples) |>
    mutate(Sex = ifelse(Sex == 1, "M", "F"))

# get sequenced samples missing from PILNK data from pedigree data and use where necessary
samples_sex <- samples_sex |>
    full_join(pedigree, by = "samples") |>
    filter(samples %in% sample_list$samples) |>
    mutate(Sex = ifelse(is.na(Sex), PedigreeSex, Sex),
           samples_bam = paste0(samples, ".bam"))

# write to samples sex file for variant calling
samples_sex |>
    dplyr::select(samples_bam, Sex) |>
    write_tsv("/mnt/loki/johnston/rumdeer/james/data/sample_sex.tsv", col_names = F)