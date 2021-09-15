ploidetect_preprocess <- function(all_data, centromeres = F, debugPlots = F, verbose = F, simplify = T, simplify_size = 500000){
  if(verbose & simplify){
    print("Ploidetect - Detection of tumour purity and aneuploidy from whole-genome sequence data")
    print("Thank you for using Ploidetect! Please remember to cite this tool if used in your research ^_^")
    print("Beginning data pre-processing steps")
  }
  # Load data
  x <- as.data.frame(all_data)
  

    
  
  # Test if data is configured and input properly
  if(any(grepl("chr", x$chr))){
    stop("Expected numeric chromosomes (1, 2, 3, ... X), not chr1, chr2, etc")
  }
  if(!all(is.numeric(x$normal))){
    stop("At least one element in normals column is not numeric.")
  }
  if(!all(is.numeric(x$tumour))){
    stop("At least one element in tumour column is not numeric.")
  }
  #if(!all(is.numeric(x$maf) | is.na(x$maf))){
  #  stop("VAF column must contain only numeric or NA values")
  #}
  if(!all(is.numeric(x$gc))){
    stop("At least one element in GC-content column is not numeric")
  }
  
  ## Step 1: Merge data such that window size is about 100Kb
  
  ## Filter for chr1-22 and X
  
  #x <- x[grepl(pattern = paste0(c(paste0("^", 1:22), "^X"), "_", collapse = "|"), x = x[,window_id]),]
  
  if(verbose){
    print("Filtering for chromosomes 1-22 and X")
  }
  x <- x %>% filter(chr %in% paste0(c(1:22, "X")))
  

  

  # Process centromere data
  if(centromeres != F){
    if(verbose){
      print("Filtering out centromeric loci")
    }
    x <- split(x, x$chr)
    centromeres_preprocess <- centromeres %>% group_by(chr) %>% dplyr::summarise(pos = first(pos), end = last(end))
    centromeres_split <- split(centromeres_preprocess, centromeres_preprocess$chr)
    x <- lapply(x, function(k){
      chr <- k$chr[1]
      centro_start <- centromeres_split[[chr]]$pos %>% unlist()
      centro_end <- centromeres_split[[chr]]$end %>% unlist()
      #print(str(k))
      k <- k %>% filter(end < centro_start | pos > centro_end)
      #print(str(k))
      return(k)
    })
    x <- rbindlist(x) %>% arrange(chr, pos)
    if(verbose){
      print("Completed centromere filtering")
    }
  }

  
  x$window_size <- x$end - x$pos
  
  
  #mean_size <- mean(x[,window_size])
  mean_size <- mean(x$window_size)
  
  closest <- round(simplify_size/mean_size, digits = 0)
  if(simplify){
    # Compute the closest integer multiple of the window size to 100kb
    if(closest < 1){
      closest = 1
    }
    # Create a merge vector
    x$merge <- floor(seq(from = 0, by = 1/closest, length.out = nrow(x)))
  }else{x$merge <- 1:nrow(x)}
  
  # Process the allele frequency column into a numeric vector
  
  #x[,avg_allele_freq][which(x[,avg_allele_freq] == ".")] <- NA
  x$maf[which(x$maf == ".")] <- NA
  
  
  #x[,avg_allele_freq] <- as.numeric(x[,avg_allele_freq])
  #x$maf <- as.numeric(x$maf)
  
  #x$chr <- gsub("_.*", "", x[,window_id])
  
  # Sanitize the column names
  
  #x <- x[,c(tumour, normal, avg_allele_freq, window_id, window_size, GC, 7, 8)]
  
  #names(x) <- c("tumour", "normal", "maf", "wind", "size", "gc", "merge", "chr")
  #x <- x %>% group_by(merge, chr) %>% summarise(tumour = sum(tumour), normal = sum(normal), maf = merge_mafs(maf, na.rm = T), wind = dplyr::first(wind), size = sum(size), gc = mean(gc))
  if(simplify){
    x <- data.table(x)
    x <- x[,.(pos = first(pos), end = last(end), tumour = sum(tumour), normal = sum(normal), maf = merge_mafs(maf, na.rm = T, exp = T), gc = mean(gc), window_size = sum(window_size)), by = list(chr, merge)]
  }


  
  ## Measure the read depth at the highest density of read coverage
  d <- density(x$tumour, n = nrow(x))
  maxpeak <- d$x[which.max(d$y)]
  
  maxpeak_ind <- which(x$tumour < (maxpeak + d$bw) & x$tumour > (maxpeak - d$bw))

  ## Remove tibble formatting
  x <- as.data.frame(x)
  
  ## Set row names to window_ids
  #row.names(x) <- x$wind

  ## Get median normal coverage
  median_normal <- median(x$normal)
  
  norm_factor <- x$normal/median_normal
  
  x$tumour <- x$tumour/norm_factor
  

  ## Perform basic pre-filtering, find the windows within the 90th percentile of tumour read counts
  #rangedf <- x[findInterval(x$tumour, 
  #                          quantile(x$tumour, 
  #                                   probs=c(0,0.90))) == 1,]
  
  ## Obtain the range of values observed within the 90th percentile of tumour read counts
  #range <- range(rangedf$tumour)[2] - range(rangedf$tumour)[1]
  
  ## Maximum allowable read depth for retention is the 90th percentile + the above computed range.
  #max <- range(rangedf$tumour)[2] + range
  #min <- 0
  
  ## Set outliers aside for later steps (these will be CNA called later, but are exempt from TC/Ploidy analysis)
  #highoutliers <- x[findInterval(x$tumour, c(min, max)) > 1,]
  
  ## Filter data for everything within the read depth range
  #x <- x[findInterval(x$tumour, c(min, max)) == 1,]
  
  # Extract X chromosome regions
  
  chrX <- x[x$chr == "X",]
  isMale=F
  if(which.min(abs((median(chrX$window_size) - median(x$window_size)) - c(0, median(x$window_size)))) == 2){
    x$normal[x$chr == "X"] <- x$normal[x$chr == "X"]/2
    isMale=T
  }
  if(verbose){
    if(isMale){
      print("Automated sex detection infers MALE sex")
    }else{
      print("Automated sex detection infers FEMALE sex")
    }
  }

  ## This is a very broad-strokes filtering step for window size. Basically removing extreme outliers w.r.t germline mappability, as we don't want to use these in modeling
  #x <- x[(x$window_size > (median(x$window_size)/10)) & (x$window_size < (median(x$window_size)*5)),]
  
  
  #x <- x[,c(tumour, normal, window_id, avg_allele_freq, window_size, GC)]
  #names(x) <- c("y_raw", "x_raw", "window", "maf", "size", "GC")
  
  #x <- x[,3:8]
  
  #names(x) <- c("y_raw", "x_raw", "maf", "window", "size", "GC")
  x <- x %>% dplyr::rename("y_raw" = "tumour", "x_raw" = "normal")
  
  plot_limits_size <- quantile(x$window_size, probs = c(0.005, 0.995))
  plot_limits_depth <- quantile(x$y_raw, probs = c(0.99))
  plot_limits_depth <- c(0, plot_limits_depth)
  if(debugPlots | !exists("debugPlots")){
    rawPlot <- x %>% ggplot(aes(x = window_size, y = y_raw)) + geom_point(size = 0.1, alpha = 0.1) + theme_bw() + xlab("Bin Size") + ylab("Tumour Read counts") + ggtitle("Raw counts by bin size")+ 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      ) + scale_x_continuous(limits = plot_limits_size) + scale_y_continuous(limits = plot_limits_depth)
    print(rawPlot)
    
    #x %>% filter(chr == 10) %>% ggplot(aes(x = pos, y = y_raw)) + geom_point(size = 0.5) + theme_bw() + ylab("Raw depth") + xlab("Position") + ggtitle("Uncorrected Read Depth")
    
    
  }
  if(debugPlots){
    GCplot <- ggplot(x, aes(x = gc * 100, y = y_raw)) + geom_point(size = 0.3, alpha = 0.2) + theme_bw() + xlab("GC Content %") + ylab("Tumour Read Counts") + ggtitle("Read count and GC content relationship") + 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      ) + scale_y_continuous(limits = plot_limits_depth) + geom_smooth(method = "lm")
    print(GCplot)  
  }
  
  ## Correct for GC content
  GCnorm <- lowesswrapper(x = x$gc, y = x$y_raw, bw = 0.75)
  if(debugPlots){
    GCnormplot <- ggplot(x, aes(y = GCnorm$residuals + maxpeak, x = window_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Corrected tumour read counts") + 
      ggtitle("Correction of read counts") + 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      ) + scale_y_continuous(limits = (plot_limits_depth - maxpeak)*1.4) + scale_x_continuous(limits = plot_limits_size)
    print(GCnormplot)
  }
  
  x$residual <- GCnorm$residuals
  
  ## Detect offset of residual value
  
  d <- density(x$residual[maxpeak_ind], n = 2^16)
  
  offset <- d$x[which.max(d$y)]
  
  x$normalized_size <- x$window_size
  
  
  
  ## Experimental scaling:
  
  #x$fan_correction <- ((x$residual/x$size) * median(x$size)) + median(x$size)
  
  
  if(debugPlots){
    GCnorm2plot <- ggplot(x, aes(y = residual, x = normalized_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Linear normalization of read counts by GC content")
    print(GCnorm2plot)
  }
  x <- x %>% mutate(residual = residual + maxpeak - offset) %>% dplyr::rename("corrected_depth" = "residual")
  d <- density(x$corrected_depth, n = 2^16)
  maxpeak <- d$x[which.max(d$y)]
  #plot(density(x$corrected_depth[x$corrected_depth < 50000])); abline(v = maxpeak)
  output <- list("x" = x, "maxpeak" = maxpeak, "merged" = closest)
  return(output)
}
