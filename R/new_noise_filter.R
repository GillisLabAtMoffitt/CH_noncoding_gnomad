# Import library
library(tidyverse)
library(GenomicRanges)

#!/usr/bin/env Rscript
# Capture the file path passed from the bash script
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]

# Example processing logic
if (!is.na(file_path)) {
  cat(paste("Currently processing:", file_path, "\n"))
  # Load your data and perform analysis here
  # data <- read.csv(file_path)
} else {
  stop("No file path provided to the script.")
}


print("args?")
args


# Load data
# gencode <- 
#   read.delim(
#     "/share/lab_gillis/Christelle/gnomAD_skweness/other_data/gencode_filtered.vcf.gz",
#     header = TRUE, sep = " ")

black_list <-
  read.delim(paste0(here::here(),
                    # "/Volumes/Lab_Gillis",
                    "/Christelle/gnomAD_skweness/analysis/data/blacklist_hg38_ENCFF356LFX.bed.gz"),
             header=FALSE) %>%
  `colnames<-`(c("seqnames", "chromStart", "chromEnd"))
# head(black_list, 2)
centromeres <-
  read.delim(paste0(here::here(),
                    # "/Volumes/Lab_Gillis",
                    "/Christelle/gnomAD_skweness/analysis/data/Centromeres_hg38.gz"))
# head(centromeres, 2)
segmental_dup <-
  read.delim(paste0(here::here(),
                    # "/Volumes/Lab_Gillis",
                    "/Christelle/gnomAD_skweness/analysis/data/Segmental_Dups_hg38.gz"))
# head(segmental_dup, 2)
wm_sdust <-
  read.delim(paste0(here::here(),
                    # "/Volumes/Lab_Gillis",
                    "/Christelle/gnomAD_skweness/analysis/data/WM_SDust_hg38.gz"))
# head(wm_sdust, 2)

gnomad <- read.delim(file_path,
                     header = TRUE
)

print(paste0("Running chr", unique(gnomad$X.CHROM)))

print("Extract AB bins")

gnomad <- gnomad %>%
  mutate(chrom = str_remove(X.CHROM, "chr"), .before = X.CHROM) %>%
  unite(IDs, c(chrom, POS, REF, ALT), sep = "-", remove = FALSE) %>%
  # 1. Extract count per AB bin for heterozygous ind--
  # ab_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for AB in heterozygous individuals; 
  # bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00>
  # mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
  separate(col = ab_hist_alt_bin_freq,
           into = c("0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
                    "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
                    "0.40-0.45","0.45-0.50","0.50-0.55","0.55-0.60",
                    "0.60-0.65","0.65-0.70","0.70-0.75","0.75-0.80",
                    "0.80-0.85","0.85-0.90","0.90-0.95","0.95-1.00"),
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
  mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
  mutate(total_allele_balance = rowSums(select(.,`0.00-0.05`:`0.95-1.00`), na.rm = TRUE)) 
# head(gnomad, 2)

removed_variants <- gnomad %>%
  select(IDs, total_allele_balance) %>% 
  filter(total_allele_balance <= 5) %>% 
  select(IDs)
print("Removed variants for total_allele_balance <= 5, N=")
print(nrow(removed_variants))
i <- unique(gnomad$chrom)
print(paste0("chr", i))
# print(paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/excluded/chr", 
#              i, "_removed_variants_total_allele_balance_inf_5.vcf.gz"))
write_delim(removed_variants,
            paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/excluded/chr",
                   i, "_removed_variants_total_allele_balance_inf_5.vcf.gz"))
rm(removed_variants)


print("Filter total_allele_balance > 5, Included N=")

gnomad <- gnomad %>%
  filter(total_allele_balance > 5) %>%
  select(IDs, everything())
print(nrow(gnomad))


print("Start fct")
# CREATE FUNCTION
noise_filter <- function(gnomad){
  
  # gnomad_vep <- gnomad %>%
  # # Extract SYMBOL from gnomAD VEP
  # mutate(ens_vep = str_match(INFO, ";vep=(.*?)$")[,2]) %>%
  # # Evaluate the number of vep in each vep string
  # mutate(number_of_vep = sapply(strsplit(ens_vep, ","), length)) %>%
  # # Separate each vep
  # separate(col = ens_vep, paste("ens_vep", 1:max(.$number_of_vep), sep="_"),
  #          sep = ",", remove = T, extra = "warn", fill = "right") %>%
  # # pivot longer VEPs
  # pivot_longer(cols = starts_with("ens_vep_"), names_to = "ENSVEP", values_to = "ens_vep") %>%
  # drop_na(ens_vep) %>% select(-ENSVEP)
  # 
  # ens_vep_var_names <- c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
  #                        "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
  #                        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
  #                        "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
  #                        "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
  #                        "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
  #                        "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
  #                        "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
  #                        "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
  #                        "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
  #                        "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
  #                        "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
  #                        "LoF_flags", "LoF_info")
  # 
  # gnomad_vep <- gnomad_vep %>%
  #   # Extract Consequence and SYMBOL
  #   separate(col = ens_vep, into = ens_vep_var_names,
  #            sep = "\\|", remove = F, extra = "warn", fill = "right") %>%
  #   select(IDs, SYMBOL) %>% 
  #   filter(!is.na(SYMBOL)) %>% 
  #   distinct(IDs, .keep_all = TRUE)

  # gnomad_vep <- gnomad_vep %>%
  #   # Add protein coding gene var
  #   left_join(., gencode %>%
  #               dplyr::rename(gencode_gene_type = gene_type),
  #             by = c("SYMBOL" = "gene_name")) %>%
  #   mutate(is_gencode_protein_coding_gene = case_when(
  #     is.na(gencode_gene_type)     ~ "No",
  #     !is.na(gencode_gene_type)    ~ "Yes"
  #   ))
  
  # # Merge with initial gnomAD
  # gnomad <- left_join(gnomad,
  #                     gnomad_vep,
  #                     by = "IDs")
  # # Cleaning
  # rm(gnomad_vep)
  
  
  
  
  # Create filter for removing noise in AB
  gnomad <- gnomad %>%
    # mutate(ab_hist_alt_bin_freq = str_match(INFO, "ab_hist_alt_bin_freq=(.*?);")[,2]) %>%
    # separate(col = ab_hist_alt_bin_freq,
    #          into = c("0.00-0.05"),
    #          sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>%
    # mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>%
    mutate(contaminated_AB = case_when(
      `0.00-0.05` == total_allele_balance               ~ "sequencing noise",
      TRUE                                              ~ "clean sequencing"
    ))# %>%
    # select(-c("0.00-0.05"))
  
  # Add cleaning lists filter
  dat <- gnomad %>% select(X.CHROM, chromStart = POS) %>% 
    mutate(chromEnd = chromStart + 0) %>% 
    mutate(seqnames = X.CHROM)
  
  dat <- makeGRangesFromDataFrame(dat)
  df <- makeGRangesFromDataFrame(black_list,
                                 start.field="chromStart",
                                 end.field="chromEnd")
  black_list <- countOverlaps(dat, df) %>% as_tibble()

  df <- makeGRangesFromDataFrame(centromeres,
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  centromeres <- countOverlaps(dat, df) %>% as_tibble()

  df <- makeGRangesFromDataFrame(segmental_dup,
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  segmental_dup <- countOverlaps(dat, df) %>% as_tibble()

  df <- makeGRangesFromDataFrame(wm_sdust,
                                 seqnames.field = "chrom",
                                 start.field="chromStart",
                                 end.field="chromEnd")
  wm_sdust <- countOverlaps(dat, df) %>% as_tibble()
  
  gnomad <- gnomad %>%
    bind_cols(., black_list %>% dplyr::rename(black_list = value),
              centromeres %>% dplyr::rename(centromeres = value),
              segmental_dup %>% dplyr::rename(segmental_duplication = value),
              wm_sdust %>% dplyr::rename(wm_sdust = value))
  
  # Add filter for depth


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  removed_variants <- gnomad %>%
    select(IDs, contaminated_AB, black_list, centromeres, segmental_duplication, wm_sdust#, dp_het_median
           ) %>% 
    filter(contaminated_AB == "sequencing noise" |
             black_list == 1 | centromeres == 1 |
             segmental_duplication != 0 |
             wm_sdust == 1 #|
             # dp_het_median < allele_depth
    )
  
  print("Removed variants for AB contamination, being in lists or low depth, N=")
  print(nrow(removed_variants))
  i <- unique(gnomad$chrom)
  print(paste0("chr", i))
  # print(paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/excluded/chr", 
  #              i, "_removed_variants_lists_or_depth.vcf.gz"))
  write_delim(removed_variants,
              paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/excluded/chr",
                     i, "_removed_variants_abcontamination_lists_or_depth.vcf.gz"))
  rm(removed_variants)
  
  # allele_depth <- 10
  gnomad <- gnomad %>%
    filter(contaminated_AB == "clean sequencing") %>%
    filter(black_list == 0) %>%
    filter(centromeres == 0) %>%
    filter(segmental_duplication == 0) %>%
    filter(wm_sdust == 0) %>% 
    # filter(dp_het_median >= allele_depth) %>% 
    select(IDs, everything())
  
}

gnomad <- noise_filter(gnomad = gnomad)
print("Filter contaminated_AB, black_list, centromeres, segmental_duplication, wm_sdust, dp_het_median")
print(paste0("Final included IDs in chr", i))
print(nrow(gnomad))
write_delim(gnomad, paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/chr",
                           i, "_noise_filter.vcf.gz"))
###################################################################### II ### Add var necessary for filtering noise variants
# for (i in chr_vec){
#   
#   ### I ### List and load files
#   file_list_gnomad <- list.files(path = "/share/lab_gillis/Christelle/gnomAD_raw_data/genome_gnomad/small_file_for_testing",
#                                  pattern = i,
#                                  recursive=FALSE,
#                                  full.names = TRUE)
#   
#   # gnomad <- do.call("rbind",lapply(Sys.glob(file_list_gnomad), read.delim,
#   #                                  header = TRUE, sep = " ",
#   #                                  colClasses=c("X.CHROM"="character"))) %>% 
#   #   mutate(X.CHROM1 = X.CHROM,
#   #          X.CHROM = str_match(X.CHROM1, "^(\\d+|X|Y)")[,1]) 
# 
#   ### II ### CALL FUNCTION
#   noise_filter(gnomad = gnomad#, cosmic = cosmic
#               )
# 
#   # n = n + 1
# 
# }


