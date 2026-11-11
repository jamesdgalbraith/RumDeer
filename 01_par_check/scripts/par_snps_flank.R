library(tidyverse)
library(plyranges)
library(Biostrings)

# read in deer seq
deer_seq <- Biostrings::readDNAStringSet("/mnt/loki/johnston/rumdeer/james/data/GCF_910594005.1_mCerEla1.1_genomic.fna")
names(deer_seq) <- sub(" .*", "", names(deer_seq))

# read in bim and pseudoautosomal snps
all_bim_snps <- read_tsv("/mnt/loki/johnston/rumdeer/james/data/20240927_Deer_1to43_mCelEla1.bim", col_names = c("chr_no", "Locus_Name", "cm", "start_20240927", "ref", "alt"))

# read in all snps and flanks
snps_flanks <- read_table("/mnt/loki/johnston/rumdeer/james/01_par_check/data/Deer_SNP_FASTA/3_DeerSequences.txt")

# read in par snps
par_snps <- read_table("/mnt/loki/johnston/rumdeer/james/data/PseudoAutosomalSNPs.txt", col_names = "Locus_Name")

# get par snps
all_snps_flanks <- snps_flanks |>
    mutate(left_flank_seq  = sub("\\[.*", "", Sequence),
           right_flank_seq = sub(".*\\]", "", Sequence))

all_snps_from_snp_data <- snps_flanks |>
    mutate(snp_data_ref = sub("\\].*", "", sub(".*\\[", "", Sequence)),
           snp_data_alt = sub(".*/", "", snp_data_ref),
           snp_data_ref = sub("/.*", "", snp_data_ref)) |>
    dplyr::select(Locus_Name, snp_data_ref, snp_data_alt)

# get flank seq as Biostings
left_flank_seq <- DNAStringSet(all_snps_flanks$left_flank_seq)
names(left_flank_seq) <- all_snps_flanks$Locus_Name
right_flank_seq <- DNAStringSet(all_snps_flanks$right_flank_seq)
names(right_flank_seq) <- all_snps_flanks$Locus_Name

# write flank seq to file
writeXStringSet(left_flank_seq, "/mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/4_DeerSNPLeftFlank_PAR.fasta")
writeXStringSet(right_flank_seq, "/mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/5_DeerSNPRightFlank_PAR.fasta")

# search genome for flank seq using blast
system("blastn -query /mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/4_DeerSNPLeftFlank_PAR.fasta -db /mnt/loki/johnston/rumdeer/james/data/GCF_910594005.1_mCerEla1.1_genomic.fna -outfmt 6 -out /mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/6_DeerSNPLeftFlank_PARBlast.tsv -num_threads 8")
message("left done")
system("blastn -query /mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/5_DeerSNPRightFlank_PAR.fasta -db /mnt/loki/johnston/rumdeer/james/data/GCF_910594005.1_mCerEla1.1_genomic.fna -outfmt 6 -out /mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/7_DeerSNPRightFlank_PARBlast.tsv -num_threads 8")
message("right done")

# read in blast hits
left_flank_blast <- read_tsv("/mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/6_DeerSNPLeftFlank_PARBlast.tsv",
                             col_names = c("qseqid", "sseqid", "pident", "length",
                                           "mismatch", "gapopen", "qstart", "qend",
                                           "sstart", "send", "evalue", "bitscore")) |>
                    dplyr::select(-c(mismatch, gapopen))
right_flank_blast <- read_tsv("/mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/7_DeerSNPRightFlank_PARBlast.tsv",
                              col_names = c("qseqid", "sseqid", "pident", "length",
                                            "mismatch", "gapopen", "qstart", "qend",
                                            "sstart", "send", "evalue", "bitscore")) |>
                    dplyr::select(-c(mismatch, gapopen))

# (if was necessary, would make sure only one copy of each flank present)

# ensure both flanks found (based on presence ion other blast)
# then work out strandness, flip coordinates where possible
left_flank_blast <- left_flank_blast |>
    dplyr::filter(qseqid %in% right_flank_blast$qseqid) |>
    mutate(left_strand = ifelse(sstart < send, "+", "-"),
           left_start =  ifelse(sstart < send, sstart, send),
           left_end   =  ifelse(sstart < send, send, sstart)) |>
    dplyr::rename(left_sseqid = sseqid) |>
    select(-c(pident, evalue, bitscore)) |>
    dplyr::rename(l_qstart = qstart, l_qend = qend, l_sstart = sstart, l_send = send, l_length = length)
right_flank_blast <- right_flank_blast |>
    dplyr::filter(qseqid %in% left_flank_blast$qseqid) |>
    mutate(right_strand = ifelse(sstart < send, "+", "-"),
           right_start =  ifelse(sstart < send, sstart, send),
           right_end   =  ifelse(sstart < send, send, sstart)) |>
    dplyr::rename(right_sseqid = sseqid) |>
    select(-c(pident, evalue, bitscore)) |>
    dplyr::rename(r_qstart = qstart, r_qend = qend, r_sstart = sstart, r_send = send, r_length = length)

# join flank blast, filter to ensure hits on same strand, calculate start, end according to strand, calulate distance between flank hits
joined_blast_hits <- inner_join(left_flank_blast, right_flank_blast, by = "qseqid") |>
    filter(left_sseqid == right_sseqid,
           left_strand == right_strand) |>
    dplyr::select(-c(right_sseqid,
                     right_strand)) |>
    dplyr::rename(seqnames = left_sseqid, strand = left_strand) |>
    dplyr::mutate(start = ifelse(strand == "+", left_end, right_end),
                  end = ifelse(strand == "+", right_start, left_start)) |>
    mutate(distance = end-start)
    # dplyr::select(qseqid, seqnames, start, end, strand, distance)

multiple_hits <- joined_blast_hits |>
                 count(qseqid) |>
                 filter(n>1)

multiple_hits_resolved <- joined_blast_hits |>
    filter(qseqid %in% multiple_hits$qseqid,
           distance == 2) |>
    dplyr::mutate(snp_pos = start + 1) |>
    dplyr::select(qseqid, seqnames, snp_pos, start, end, strand)

single_hits_resolved <- joined_blast_hits |>
    filter(!qseqid %in% multiple_hits$qseqid) |>
    filter(distance == 2) |>
    dplyr::mutate(snp_pos = start + 1) |>
    dplyr::select(qseqid, seqnames, snp_pos, start, end, strand)

# if either flank is 150bp can calculate snp_pos
resolvable_joined_blast_single_hits <- joined_blast_hits |>
    filter(!qseqid %in% multiple_hits$qseqid) |>
    filter(distance > 2) |>
    filter(l_length == 150 | r_length == 150) |>
    mutate(snp_pos = ifelse(l_length == 150, start + 1, end - 1)) |>
    dplyr::select(qseqid, seqnames, snp_pos, start, end, strand)

# resolved_hits 
compiled_resolved_hits <- single_hits_resolved |>
    rbind(multiple_hits_resolved) |>
    rbind(resolvable_joined_blast_single_hits) |>
    dplyr::rename(Locus_Name = qseqid)

compiled_resolved_hits |>
    inner_join(par_bim_snps) |>
    select(-c(strand, chr_no, cm, snp_pos, alt))

# 
unresolved_hits <- snps_flanks[!snps_flanks$Locus_Name %in% compiled_resolved_hits$Locus_Name,]

# get resolved hits in PAR
resolved_PAR_snps <- compiled_resolved_hits |>
    filter(Locus_Name %in% par_snps$Locus_Name)

# get resolved hits in PAR
resolved_PAR_snp_ranges <- resolved_PAR_snps |>
    dplyr::mutate(start = snp_pos) |>
    mutate(width = 1) |>
    plyranges::as_granges()
resolved_PAR_snp_seq <- BSgenome::getSeq(deer_seq, resolved_PAR_snp_ranges)
names(resolved_PAR_snp_seq) <- resolved_PAR_snp_ranges$Locus_Name

# combined all data together (i.e. data fom PLINK .bim, 3_DeerSequences.txt and the sequence found using BLAST)
all_data_combined <- all_bim_snps |>
    inner_join(
        all_snps_from_snp_data,
        by = "Locus_Name") |>
    inner_join(
        tibble(actual_seq = as.character(resolved_PAR_snp_seq), Locus_Name = names(resolved_PAR_snp_seq)),
        by = "Locus_Name") |>
    inner_join(
        resolved_PAR_snps) |>
    dplyr::select(-c(cm, chr_no))
all_data_combined |>
    write_tsv("/mnt/loki/johnston/rumdeer/james/01_par_check/outputs/PAR_search/combined_output.tsv")