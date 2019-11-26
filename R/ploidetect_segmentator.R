ploidetect_segmentator <- function(processed_data, maxpeak, verbose = F, segmentation_threshold = segmentation_threshold, tp = tp, ploidy = ploidy){
  if(verbose){
    print("Performing segmentation and copy number variant calling")
  }
  # Start with all NA values for CN
  processed_data$CN <- NA
  # Extract relevant data for CNA calling
  data <- processed_data[,c("chr", "pos", "end", "maf", "corrected_depth", "CN")]
  # Extract windows that are extreme outliers in coverage and include them in CNA calling (since they were filtered out in data pre-processing for TC/ploidy modeling)
  ## Split data by chromosome
  datsplit <- split(data, data$chr)
  
  ## Generate predicted positions of each CN state given the TP and Ploidy

  cov_char <- get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)
  
  predictedpositions <- cov_char$cn_by_depth
  predicted_homd <- cov_char$homd
  predicted_diff <- cov_char$diff
  
  
  ## Generate model for regression-based CNA calling
  df_train <- data.frame("CN" = as.numeric(names(predictedpositions)), "segment_depth" = predictedpositions, stringsAsFactors = F)
  train <- lm(CN ~ segment_depth, data = df_train)
  ## Run segmentation by compression on all chromosomes
  if(verbose){
    print("Performing segmentation of copy number data")
  }
  compressedalldat <- lapply(datsplit, runiterativecompression, x = predicted_diff, segmentation_threshold = segmentation_threshold, verbose = verbose)
  
  ## Generate a "median_segment" column for median coverage per segment
  compressedalldat <- lapply(compressedalldat, function(x){
    #start <- Sys.time()
    segs <- x$segment
    x <- data.table(x)
    x$segment_depth <- x[,median(corrected_depth), by = x$segment][segs][,2]
    x$median_maf = x[,merge_mafs(maf, exp = T), by = segs][segs][,2]
    return(x)
    #x$median_maf = unlist(tapply(x$maf, x$segment, merge_mafs, na.rm = T, exp = T))
    #end <- Sys.time()
    #print(end-start)
    #start <- Sys.time()
    #x <- x %>% group_by(segment) %>% dplyr::mutate("segment_depth" = median(corrected_depth), "median_maf" = merge_mafs(maf, na.rm = T, exp = T))
    #end <- Sys.time()
    #print(end-start)
  })
  return(compressedalldat)
}
