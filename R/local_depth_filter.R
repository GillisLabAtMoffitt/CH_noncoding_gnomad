## Import library

library(tidyverse)


gnomad <- read.delim(paste0(here::here(), 
                            "/raw_data/chr2_noise_filter_DPbinExtract_limdata.txt"), 
                     header = TRUE)
gnomad <- gnomad %>% 
  # slice(1:10) |> 
  mutate(chrom = str_remove(X.CHROM, "chr"), .before = X.CHROM) %>%
  unite(IDs, c(chrom, POS, REF, ALT), sep = "-", remove = FALSE)

print(paste0("Running", unique(gnomad$X.CHROM)))


# CREATE FUNCTION
# depth_extract <- function(gnomad){

  gnomad_dp <- gnomad %>% 
    
    # select(IDs, INFO) %>% 
    mutate(IDs = as.character(IDs)) %>% 
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
  
  
  gnomad_dp <- gnomad_dp %>% 
    pivot_longer(cols = c(
      starts_with("dp_")
    ), 
    names_pattern = "(.*)_(...)$", 
    names_to = c("depth_bin", "depth_type"),
    values_to = "dp_value") %>% 
    mutate(depth_bin = str_remove(depth_bin, "dp_"))
  
  
  
  gnomad_dp <- gnomad_dp %>%
    mutate(depth_bin_center = case_when(
      depth_bin == "0-5"              ~ 2.5,
      depth_bin == "5-10"             ~ 7.5,
      depth_bin == "10-15"            ~ 12.5,
      depth_bin == "15-20"            ~ 17.5,
      depth_bin == "20-25"            ~ 22.5,
      depth_bin == "25-30"            ~ 27.5,
      depth_bin == "30-35"            ~ 32.5,
      depth_bin == "35-40"            ~ 37.5,
      depth_bin == "40-45"            ~ 42.5,
      depth_bin == "45-50"            ~ 47.5,
      depth_bin == "50-55"            ~ 52.5,
      depth_bin == "55-60"            ~ 57.5,
      depth_bin == "60-65"            ~ 62.5,
      depth_bin == "65-70"            ~ 67.5,
      depth_bin == "70-75"            ~ 72.5,
      depth_bin == "75-80"            ~ 77.5,
      depth_bin == "80-85"            ~ 82.5,
      depth_bin == "85-90"            ~ 87.5,
      depth_bin == "90-95"            ~ 92.5,
      depth_bin == "95-100"           ~ 97.5,
      TRUE                            ~ 105
    )) %>% 
    mutate(dp_frequency = case_when(
      depth_type == "all"               ~ dp_value / nbr_individuals_with_dp,
      depth_type == "het"               ~ dp_value / nbr_het_individuals_with_dp
    ))

  dp_all_median <-  data.frame(matrix(nrow=length(unique(gnomad_dp$IDs)), ncol=1,
                                      dimnames = list(c(unique(gnomad_dp$IDs)),
                                                      c("value")) ))
  dp_het_median <-  data.frame(matrix(nrow=length(unique(gnomad_dp$IDs)), ncol=1,
                                      dimnames = list(c(unique(gnomad_dp$IDs)),
                                                      c("value")) ))
  tictoc::tic
  set.seed(1234)
  for(i in unique(gnomad_dp$IDs)) {
    
    data_plot <- gnomad_dp %>% 
      filter(IDs == i
      )
    
    p <- data_plot %>% 
      ggplot(aes(x = depth_bin_center, y= dp_frequency, color= depth_type
      ))+
      geom_smooth(alpha= 0.5, se = FALSE, method = "loess",
                  span = 0.6
      )+
      ylim(0, max(data_plot$dp_frequency))
    p_ <- layer_data(p, 1) %>% 
      mutate(colour =  as.factor(colour)) %>% 
      mutate_if(is.numeric, ~replace(., is.na(.), 0))
    
    all <- p_ %>% filter(group == 1) 
    all <- sample(all$x, 180, prob = all$y, replace = TRUE)
    variant <- p_ %>% filter(group == 2) 
    variant <- sample(variant$x, 180, prob = variant$y, replace = TRUE)
    
    dp_all_median[i,] <- median(all, na.rm = TRUE)
    dp_het_median[i,] <- median(variant, na.rm = TRUE)

  }
  tictoc::toc()
  rm(gnomad_dp)

  dp_all_median <- dp_all_median %>% 
    bind_cols(., dp_het_median) %>% 
    `colnames<-`(c("dp_all_median", "dp_het_median")) %>% 
    rownames_to_column("IDs")
  
  # Bind the distribution characteristic values
  gnomad <- gnomad %>%
    full_join(., dp_all_median, by = "IDs")
  rm(dp_all_median, dp_het_median)
  
  removed_variants <- gnomad %>% 
    filter(dp_het_median < 10)
i <- unique(gnomad$chrom)
  write_delim(removed_variants,
              paste0(here::here(), "/data/processed_data",
                     "/chr",
                     i, "_removed_variants_abcontamination_lists_or_lowdepth.txt"))
  
  gnomad <- gnomad %>% 
    filter(dp_het_median >= 10)
  
# }


########## Run
# gnomad <- depth_extract(gnomad = gnomad)
print(paste0("Final included IDs in chr", i))
print(nrow(gnomad))
write_delim(gnomad, 
            paste0(here::here(), "/data/processed_data",
                   "/chr",
                   i, "_depth_filter.txt"))

# for (i in chr_vec){
#   
#   ### I ### List and load files
#   
#   file_list <- list.files(path = "/share/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/noise_filter_full_data",
#                           # path = "/Volumes/lab_gillis/Christelle/intronic_skewness/pipeline/noise_filter/noise_filter_full_data",
#                           pattern = i,
#                           recursive = FALSE,
#                           full.names = TRUE)
#   
#   ### II ### Run fct
#   depth_extract(gnomad = gnomad)
#   
# 
# }




