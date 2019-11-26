ploidetect_fineCNAs <- function(all_data, CNAout, TC, ploidy, depthdiff = depthdiff, maxpeak = maxpeak, verbose = verbose, decision = decision, simpsize = simpsize, unaltered = unaltered){
  ## Generate processed data.frame for unmerged data
  unmerged_data <- ploidetect_preprocess(all_data = all_data, verbose = verbose, debugPlots = F, simplify = T, simplify_size = simpsize/2)

  den <- density(unmerged_data$x$corrected_depth, n = nrow(unmerged_data$x))
  offset <- den$x[which.max(den$y)]
  unmerged_data$x$corrected_depth <- unmerged_data$x$corrected_depth - offset
  
  ## Unpack maxpeak
  unmerged_maxpeak <- unmerged_data$maxpeak
  
  ## Compute reads-per-copy and HOMD location for unmerged data based on merged data
  unmerged_diff <- depthdiff/(unaltered/unmerged_data$merged)
  unmerged_normalreads <- unmerged_maxpeak - ploidy*unmerged_diff
  
  predictedpositions <- seq(from = unmerged_normalreads, by = unmerged_diff, length.out = 11)
  names(predictedpositions) <- 0:10
  
  df.train <- data.frame("CN" = 0:10, "median_segment" = predictedpositions)
  model <- lm(CN ~ median_segment, data = df.train)
  #print(predictedpositions)
  
  ## Continue unpacking data
  #unmerged_highoutliers <- unmerged_data$highoutliers %>% dplyr::rename("y_raw" = "tumour", "x_raw" = "normal") %>% mutate("residual" = y_raw, "normalized_size" = window_size)
  unmerged_data <- unmerged_data$x
  if(decision == 2){
    unmerged_data$residual <- unmerged_data$y_raw - unmerged_maxpeak
  }
  
  ## Sanitise highoutliers and merge with rest of data
  #unmerged_highoutliers <- unmerged_highoutliers[,c("tumour", "normal", "maf", "wind", "size", "gc", "tumour", "size")]
  #names(unmerged_highoutliers) <- names(unmerged_data$x)
  unmerged_data$corrected_depth <- unmerged_data$corrected_depth + unmerged_maxpeak
  #unmerged_data <- rbind.data.frame(unmerged_data, unmerged_highoutliers)
  
  ## Compute the positional columns (chr, pos, end) for each window
  #unmerged_data$chr <- gsub("_.*", "", unmerged_data$window)
  #unmerged_data$pos <- as.numeric(gsub(".*_", "", unmerged_data$window))
  #unmerged_data$end <- unmerged_data$pos + unmerged_data$size
  unmerged_data <- unmerged_data %>% arrange(chr, pos)
  
  

  
  
  ## Split both merged and unmerged data by chr
  unmerged_data <- split(unmerged_data, f = unmerged_data$chr)
  merged_data <- split(CNAout, f = CNAout$chr)
  
  ## First map old CNs to new higher-res data
  for(chr in names(unmerged_data)){
    unmerged_chr <- unmerged_data[[chr]] %>% arrange(pos)
    merged_chr <- merged_data[[chr]]
    merged_segments <- merged_chr %>% group_by(chr, segment) %>% dplyr::summarise(pos = dplyr::first(pos), end = last(end), CN = mean(CN), maf = merge_mafs(maf, exp = T)) %>% arrange(pos)
    merged_segments$pos[1] <- 0
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, merged_segments$pos)
    unmerged_chr <- left_join(unmerged_chr, merged_segments[,c("segment", "CN")], by = "segment")
    unmerged_data[[chr]] <- unmerged_chr
  }
  
  #unmerged_data$`11` %>% ggplot(aes(x = pos,  y = residual, color = segment)) + geom_point() + scale_color_viridis()

  ## Compute the standard deviation of read depth in the 50% longest segments
  grouped_data <- do.call(rbind.data.frame, unmerged_data)
  sd <- grouped_data %>% group_by(chr, segment) %>% dplyr::summarise("sd" = sd(corrected_depth), "mean_residual" = mean(corrected_depth), "length" = n()) %>% ungroup %>% arrange(desc(length)) %>% slice(n()/2) %>% summarise("medsd" = median(sd, na.rm = T)) %>% unlist
  #print(sd)
  #test_data$new_CN <- round(predict(model, data.frame("median_segment" = test_data$residual)), 0)
  #test_data$flagged <- test_data$CN != test_data$new_CN
  #test_data %>% filter(chr == 1, residual < 1e+06) %>% ggplot(aes(x = pos, y = residual, color = flagged)) + geom_point() + scale_color_viridis(discrete = T)
  
  #unmerged_data$`1` %>% ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis()
  
  for(chr in names(unmerged_data)){
    unmerged_chr <- unmerged_data[[chr]] %>% arrange(pos)
    merged_chr <- merged_data[[chr]]
    merged_segments <- merged_chr %>% group_by(chr, segment) %>% summarise(pos = dplyr::first(pos), end = last(end), CN = mean(CN))
    merged_segments$pos[1] <- 0
    #unmerged_chr$segment <- findInterval(unmerged_chr$pos, merged_segments$pos)
    #unmerged_chr <- left_join(unmerged_chr, merged_segments[,c("segment", "CN")], by = "segment")
    #unmerged_chr$mafflipped <- abs(unmerged_chr$maf - 0.5) + 0.5
    #unmerged_chr$new_CN <- round(predict(model, data.frame("median_segment" = unmerged_chr$residual)), 0)
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% dplyr::mutate("mean_residual" = mean(corrected_depth), "z" = (corrected_depth - mean(corrected_depth))/sd)
    #unmerged_chr %>% ggplot(aes(x = pos, y = corrected_depth, color = abs(z) > 3)) + geom_point() + scale_color_viridis(discrete = T)
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% dplyr::mutate("median_segment" = median(corrected_depth))
    unmerged_chr <- unmerged_chr[,c("chr", "pos", "end", "corrected_depth", "CN", "segment", "median_segment", "z", "maf")] %>% arrange(pos)
    unmerged_chr$flagged <- F
    #unmerged_chr$flagged[which(abs(unmerged_chr$residual - unmerged_chr$median_segment) > unmerged_diff * 0.5)] <- T
    unmerged_chr$flagged[which(abs(unmerged_chr$z) > 3)] <- T
    #unmerged_chr %>% ggplot(aes(x = pos, y = residual, color = flagged)) + geom_point() + scale_color_viridis(discrete = T)
    whichflagged <- which(unmerged_chr$flagged)
    for(flagged in whichflagged){
      if(!any(c(flagged - 1, flagged + 1) %in% whichflagged)){
        unmerged_chr$flagged[flagged] <- F
      }
    }
    #unmerged_chr %>% filter(CN < 10) %>%  ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    #unmerged_chr$flagged <- unmerged_chr$CN != unmerged_chr$new_CN
    unmerged_chr$old_segment <- F
    #merged_chr[which(merged_chr$breakpoint),]
    
    ## Identify intervals of old breakpoints
    #breakpoints <- c(merged_chr$pos[which(merged_chr$breakpoint)], merged_chr$end[which(merged_chr$breakpoint)]) %>% sort()
    
    breakpoints <- merged_chr %>% group_by(segment) %>% dplyr::summarise("5-term_pos" = first(pos) - 1, "5-term_end" = first(end) + 1, "3-term_pos" = last(pos) - 1, "3-term_end" = last(end) + 1)
    
    breakpoints <- as.matrix(breakpoints[,2:5]) %>% t() %>% as.vector() %>% sort()
    
    ## All windows that fall within old breakpoints need to be broken into length=1 segments
    ## First we generate intervals based on the old breakpoints
    unmerged_chr$atomize <- findInterval(unmerged_chr$pos, breakpoints)
    intervals <- unique(unmerged_chr$atomize)
    ## Find odd-numbered intervals, which denote the points which actually fell within the range of old breakpoints
    relevant_intervals <- intervals[which(intervals %% 2 == 1)]
    
    breakpoints <- c(unmerged_chr$pos[which(unmerged_chr$atomize %in% relevant_intervals)], unmerged_chr$end[which(unmerged_chr$atomize %in% relevant_intervals)]) %>% sort()
    
    #breakpoints <- c(breakpoints, breakpoints + 1)
    new_segment_interval <- unmerged_chr$pos[which(unmerged_chr$flagged)]

    new_segment_interval <- sort(c(new_segment_interval, new_segment_interval + 1))
    new_segment_interval <- sort(c(breakpoints, new_segment_interval))
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, new_segment_interval, rightmost.closed = F) + 1
    #unmerged_chr %>% filter(CN < 10) %>% ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    unmerged_to_compress <- unmerged_chr %>% dplyr::mutate("npoints" = 1) %>% group_by(segment) %>% arrange(pos) %>% dplyr::summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "corrected_depth" = sum(corrected_depth), "maf" = merge_mafs(maf, exp = T)) %>% arrange(pos) %>% dplyr::mutate("mean_residual" = corrected_depth/npoints)
    #unmerged_to_compress %>% ggplot(aes(x = pos, xend = end, y = residual/npoints, yend = residual/npoints)) + geom_segment() + geom_point(data = unmerged_chr_filt, mapping = aes(x = pos, y = residual, color = segment), inherit.aes = F) + scale_color_viridis()
    
    
    
    new_segments <- runiterativecompression(t = unmerged_to_compress, x = unmerged_diff, segmentation_threshold = 0.5, verbose = verbose)
    #new_segments <- compressdata(t = unmerged_to_compress, x = unmerged_diff, segmentation_threshold = 0.25) %>% mutate("segment" = 1:n())
    #new_segments %>% filter(CN < 10) %>% ggplot(aes(x = pos, y = residual/npoints, color = segment)) + geom_point() + scale_color_viridis() # + geom_point(data = unmerged_chr_filt, mapping = aes(x = pos, y = residual, color = segment), inherit.aes = F) + scale_color_viridis()
    
    new_segments <- new_segments %>% group_by(segment) %>% dplyr::summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "corrected_depth" = sum(corrected_depth), "maf" = merge_mafs(maf, exp = T), "len" = n())

    #long_segment <- 0
    #new_chr = F
    #if(nrow(new_segments) > 1){
    #  for(segment in 1:nrow(new_segments)){
    #    ## If there's a new chromosome, reset long_segment
    #    if(segment > 1){
    #      if(new_segments$chr[segment] != new_segments$chr[segment-1]){
    #        long_segment <- 0
    #      }
    #    }
    #    ### Lookahead to get new long_segment
    #    ## If long_segment is zero
    #    if(long_segment == 0){
    #      ## Record current chromosome
    #      current_chr = new_segments$chr[segment]
    #      ## Lookahead
    #      for(sub_segment in segment:nrow(new_segments)){
    #        ## If we get to a new chromosome, break
    #        if(current_chr != new_segments$chr[sub_segment]){
    #          new_chr = T
    #          break
    #        }
    #        ## When long_segment is found, record it
    #        if(new_segments$npoints[sub_segment] > 1){
    #          long_segment <- new_segments$segment[sub_segment]
    #          break
    #        }
    #      }
    #      if(new_chr){
    #        new_chr = F
    #        next
    #      }
    #    }
    #    ## If the current segment is "long" ie. supported by at least two adjacent points, make this the new "long_segment"
    #    if(new_segments$npoints[segment] > 1){
    #      long_segment <- new_segments$segment[segment]
    #    }
    #    ## If we haven't found a new long_segment yet, skip#

        ## If the current segment is length one, merge to long_segment
    #    if(new_segments$npoints[segment] == 1){
    #      new_segments$segment[segment] <- long_segment
    #    }
    #  }
    #}
    #new_segments <- new_segments %>% ungroup %>% group_by(segment) %>% summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "residual" = sum(residual), "len" = n())
    #print(new_segments)
    print("Generating segments based on compression data")
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, new_segments$pos, rightmost.closed = F)
    
    centromere_start <- centromeres$pos[centromeres$chr == chr][1]
    centromere_end <- centromeres$end[centromeres$chr == chr][2]
    
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% dplyr::mutate("median_segment" = median(corrected_depth))
    print("calling copy number")
    unmerged_chr$CN <- round(predict(model, unmerged_chr), digits = 1)
    print("calling breakpoints")
    unmerged_chr <- callbreakpoints(unmerged_chr, predictedpositions = predictedpositions, maxpeak = unmerged_maxpeak)
    unmerged_chr %>% filter(CN < 10) %>%  ggplot(aes(x = pos, y = corrected_depth, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    names(unmerged_chr)[which(names(unmerged_chr) == "residual")] <- "raw_residual"
    names(unmerged_chr)[which(names(unmerged_chr) == "residual")] <- "residual"
    unmerged_data[chr] <- list(unmerged_chr)
  }
  return(unmerged_data)
}


ploidetect_resegment <- function(CNAout, TC, ploidy, depthdiff = depthdiff, maxpeak = maxpeak, verbose = verbose, decision = decision){
  ## Convert input data into format for iterative compression
  merged_data <- CNAout %>% group_by(chr, segment) %>% dplyr::summarise(pos = first(pos), end = last(end), "corrected_depth" = sum(corrected_depth), "npoints" = n())
  ## Split by chromosome
  list_data <- split(merged_data, f = merged_data$chr)
  ## Run iterative compression segmentation
  compressed <- lapply(list_data, runiterativecompression, x = depthdiff, segmentation_threshold = 0.5)
  ## Plot for debugging, ignore
  #bk[[1]] %>% ggplot(aes(x = pos, y = corrected_depth/npoints, color = segment)) + geom_point() + scale_color_viridis()
  ## Turn back into dataframe
  compressed <- do.call(rbind.data.frame, compressed)
  ## Summarize by segment
  compressed <- compressed %>% group_by(chr, segment) %>% dplyr::summarise("pos" = first(pos), "end" = last(end), "n" = n())
  #compressed$pos <- unlist(tapply(compressed, list(compressed$chr, compressed$segment)), function(x)x[1,])
  ## Split by chromosome for mapping back to input data
  compressed <- split(compressed, f = compressed$chr)
  out_data <- split(CNAout, f = CNAout$chr)
  ## Map segments by chromosome
  for(i in names(out_data)){
    out_data[[i]]$segment <- findInterval(x = out_data[[i]]$pos, vec = compressed[[i]]$pos)
  }
  ## Convert back to df
  out_data <- do.call(rbind.data.frame, out_data)
  ## Compute segment statistics again
    
  out_data$segment_depth = unlist(tapply(out_data$corrected_depth, list(out_data$chr, out_data$segment), function(x)rep(median(x), times = length(x))))
  out_data$median_maf = unlist(tapply(out_data$maf, list(out_data$chr, out_data$segment), function(x)rep(merge_mafs(x, na.rm = T, exp = T), times = length(x))))
  
  ## Create training dataframe for copy number calling
  df.train <- data.frame("CN" = seq(from = ploidy - 5, to = ploidy + 5, by = 1), "segment_depth" = seq(from = maxpeak - 5*depthdiff, to = maxpeak + 5*depthdiff, length.out = 11))
  model <- lm(CN ~ segment_depth, data = df.train)
  out_data$CN <- round(predict(model, out_data), 2)
  out_data <- callbreakpoints(out_data, predictedpositions = seq(from = maxpeak - 5*depthdiff, to = maxpeak + 5*depthdiff, length.out = 11), maxpeak = maxpeak)
  return(out_data)
}
