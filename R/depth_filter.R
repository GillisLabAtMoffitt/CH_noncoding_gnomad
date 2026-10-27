## Import library

library(tidyverse)

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

# --- runtime / memory reporting (interactive HPC observability) ----------
.t0 <- proc.time()
invisible(gc(reset = TRUE, verbose = FALSE))

gnomad <- read_tsv(file_path, show_col_types = FALSE)

print(paste0("Running", unique(gnomad$X.CHROM)))

print("Extract DP bins")

gnomad <- gnomad %>% 
  mutate(chrom = str_remove(X.CHROM, "chr"), .before = X.CHROM) %>%
  unite(IDs, c(chrom, POS, REF, ALT), sep = "-", remove = FALSE) %>% 
  
  # select(IDs, INFO) %>% 
  # mutate(IDs = as.character(IDs)) %>% 
  distinct(IDs, .keep_all = TRUE) %>%
  # dp_hist_all_bin_freq,Number=A,Type=String,Description="Histogram for DP; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
  # mutate(dp_hist_all_bin_freq = str_match(INFO, "dp_hist_all_bin_freq=(.*?);")[,2]) %>%
  
  separate(col = dp_hist_all_bin_freq,
           into = c("dp_0-5_all","dp_5-10_all","dp_10-15_all","dp_15-20_all",
                    "dp_20-25_all","dp_25-30_all","dp_30-35_all","dp_35-40_all",
                    "dp_40-45_all","dp_45-50_all","dp_50-55_all","dp_55-60_all",
                    "dp_60-65_all","dp_65-70_all","dp_70-75_all","dp_75-80_all",
                    "dp_80-85_all","dp_85-90_all","dp_90-95_all","dp_95-100_all"), 
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>% 
  # dp_hist_all_n_larger,Number=A,Type=Integer,Description="Count of DP values falling above highest histogram bin edge">
  # mutate("dp_>100_all" = as.numeric(str_match(INFO, "dp_hist_all_n_larger=(.*?);")[,2])) %>% 
  rename("dp_>100_all" = dp_hist_all_n_larger) %>% 
  # mutate(across(grep("^[[:digit:]]", colnames(.)), ~ as.numeric(.))) %>% 
  
  # dp_hist_alt_bin_freq,Number=A,Type=String,Description="Histogram for DP in heterozygous individuals; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100">
  # mutate(dp_hist_alt_bin_freq = str_match(INFO, "dp_hist_alt_bin_freq=(.*?);")[,2]) %>%
  separate(col = dp_hist_alt_bin_freq,
           into = c("dp_0-5_het","dp_5-10_het","dp_10-15_het","dp_15-20_het",
                    "dp_20-25_het","dp_25-30_het","dp_30-35_het","dp_35-40_het",
                    "dp_40-45_het","dp_45-50_het","dp_50-55_het","dp_55-60_het",
                    "dp_60-65_het","dp_65-70_het","dp_70-75_het","dp_75-80_het",
                    "dp_80-85_het","dp_85-90_het","dp_90-95_het","dp_95-100_het"), 
           sep = "\\|", remove = TRUE, extra = "warn", fill = "right") %>% 
  # dp_hist_alt_n_larger,Number=A,Type=Integer,Description="Count of DP values falling above highest histogram bin edge">
  # mutate("dp_>100_het" = as.numeric(str_match(INFO, "dp_hist_alt_n_larger=(.*?);")[,2])) %>% 
  rename("dp_>100_het" = dp_hist_alt_n_larger) %>% 
  # select(-INFO) %>% 
  mutate(across(starts_with("dp_"), ~ as.numeric(.))) %>% 
  
  mutate(nbr_individuals_with_dp = rowSums(select(.,`dp_0-5_all`:`dp_>100_all`), na.rm = TRUE)) %>%
  mutate(nbr_het_individuals_with_dp = rowSums(select(.,`dp_0-5_het`:`dp_>100_het`), na.rm = TRUE))

removed_variants <- gnomad %>%
  select(IDs, nbr_het_individuals_with_dp) %>% 
  filter(nbr_het_individuals_with_dp == 0) %>% 
  select(IDs)
print("Removed variants for nbr_het_individuals_with_dp == 0, N=")
print(nrow(removed_variants))
i <- unique(gnomad$chrom)
print(paste0("chr", i))
# print(paste0("/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/excluded/chr", 
#              i, "_removed_variants_total_allele_balance_inf_5.vcf.gz"))
write_delim(removed_variants,
            paste0("excluded/chr",
                   i, "_removed_variants_nbr_het_individuals_with_dp_equal_0.vcf.gz"))
rm(removed_variants)


print("Filter nbr_het_individuals_with_dp > 0, Included N=")

gnomad <- gnomad %>%
  filter(nbr_het_individuals_with_dp > 0) %>%
  select(IDs, everything())
print(nrow(gnomad))


print("Run fct")

# Vectorised weighted median of the DP histogram bins (linear interpolation
# within the bin where the cumulative weight first crosses half). Replaces the
# original per-variant ggplot+LOESS+sample(180) loop. See CHANGELOG 2026-06-12
# for the validation against chr1/chr2 golden masters and the speedup numbers.
weighted_median_rows <- function(weights) {
  n <- nrow(weights); k <- ncol(weights)
  bin_low   <- c(seq(0, 95, by = 5), 100)[seq_len(k)]
  bin_width <- c(rep(5, 20), 10)[seq_len(k)]   # >100 bin given width 10
  totals    <- rowSums(weights, na.rm = TRUE)
  half      <- totals / 2
  utri      <- upper.tri(matrix(0, k, k), diag = TRUE) * 1
  cs        <- weights %*% utri
  cs_prev   <- cbind(0, cs[, -k, drop = FALSE])
  idx       <- max.col(cs >= half, ties.method = "first")
  rc        <- cbind(seq_len(n), idx)
  med       <- bin_low[idx] + (half - cs_prev[rc]) / weights[rc] * bin_width[idx]
  med[totals == 0] <- NA_real_
  med
}

# CREATE FUNCTION
depth_extract <- function(gnomad){
  bin_cols <- function(suffix) c(
    paste0("dp_", c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40",
                    "40-45","45-50","50-55","55-60","60-65","65-70","70-75","75-80",
                    "80-85","85-90","90-95","95-100"), "_", suffix),
    paste0("dp_>100_", suffix)
  )
  gnomad$dp_all_median <- weighted_median_rows(as.matrix(gnomad[, bin_cols("all")]))
  gnomad$dp_het_median <- weighted_median_rows(as.matrix(gnomad[, bin_cols("het")]))

  removed_variants <- gnomad %>%
    select(IDs, dp_het_median) %>%
    filter(dp_het_median < 10)
  print("Removed variants for dp_het_median < 10, N=")
  print(nrow(removed_variants))
  i <- unique(gnomad$chrom)
  write_delim(removed_variants,
              paste0("excluded/chr",
                     i, "_removed_variants_depth_filter.vcf.gz"))

  gnomad <- gnomad %>%
    filter(dp_het_median >= 10) %>%
    select(IDs, chrom, X.CHROM, POS, REF, ALT,
           dp_all_median, nbr_individuals_with_dp, dp_het_median, nbr_het_individuals_with_dp) %>%
    distinct(IDs, .keep_all = TRUE)

}


########## Run
gnomad <- depth_extract(gnomad = gnomad)
i <- unique(gnomad$chrom)
print(paste0("Final included IDs in chr", i))
print(nrow(gnomad))
write_delim(gnomad, paste0("chr",
                           i, "_depth_filter.txt"))

# --- end-of-run report ----------------------------------------------------
.dt <- proc.time() - .t0
.heap_mb <- sum(gc(verbose = FALSE)[, ncol(gc(verbose = FALSE))])
.rss_mb <- tryCatch(suppressWarnings({
  hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(sub(".*?([0-9]+).*", "\\1", hwm)) / 1024
}), error = function(e) NA_real_)
cat(sprintf(
  "\n=== runtime: %.1fs (user=%.1f sys=%.1f) | R heap peak: %.0f MB | process RSS peak: %s ===\n",
  .dt[["elapsed"]], .dt[["user.self"]], .dt[["sys.self"]], .heap_mb,
  if (is.na(.rss_mb)) "n/a (non-Linux)" else sprintf("%.0f MB", .rss_mb)))




