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
  if(!all(is.numeric(x$gc))){
    stop("At least one element in GC-content column is not numeric")
  }
  
  ## Step 1: Merge data such that window size is about 100Kb
  
  ## Filter for chr1-22 and X

  
  if(verbose){
    print("Filtering for chromosomes 1-22 and X")
  }
  x <- x %>% filter(chr %in% paste0(c(1:22, "X")))
  

  

  # Process centromere data
  if(centromeres != F){
    if(verbose){
      print("Filtering out centromeric loci")
    }
    ## Split by chromosome
    x <- split(x, x$chr)
    ## Strip "chr" from chromosomes if present
    if(grepl("chr", centromeres$chr[1])){
      centromeres[,chr:=gsub("chr", "", chr)]
    }
    ## Not entirely sure why this is here, 
    ## TODO: Test this
    centromeres_preprocess <- centromeres %>% group_by(chr) %>% dplyr::summarise(pos = first(pos), end = last(end))
    ## Split centromeres by chromosome
    centromeres_split <- split(centromeres_preprocess, centromeres_preprocess$chr)
    ## Loop by chromosome, filter out centromeric positions
    x <- lapply(x, function(k){
      chr <- k$chr[1]
      centro_start <- unlist(centromeres_split[[chr]]$pos)
      centro_end <- unlist(centromeres_split[[chr]]$end)
      k <- k %>% filter(end < centro_start | pos > centro_end)
      return(k)
    })
    ## Merge data back from list form
    x <- rbindlist(x) %>% arrange(chr, pos)
    if(verbose){
      print("Completed centromere filtering")
    }
  }

  ## Get mean bin size
  x$window_size <- x$end - x$pos
  mean_size <- mean(x$window_size)
  ## Figure out how many times we need to merge bins to create low-resolution TC/Ploidy bins
  closest <- round(simplify_size/mean_size, digits = 0)
  if(simplify){
    # Compute the closest integer multiple of the window size to the simplify size
    if(closest < 1){
      closest = 1
    }
    # Create a merge vector
    x$merge <- floor(seq(from = 0, by = 1/closest, length.out = nrow(x)))
  }else{x$merge <- 1:nrow(x)}
  
  # Convert "." in baf column to NA
  x$maf[which(x$maf == ".")] <- NA
  ## Merge based on the aggregation column
  if(simplify){
    x <- data.table(x)
    x <- x[,.(pos = first(pos), end = last(end), tumour = sum(tumour), normal = sum(normal), maf = merge_mafs(maf, na.rm = T, exp = T), gc = mean(gc), window_size = sum(window_size)), by = list(chr, merge)]
  }


  
  ## Measure the read depth at the highest density of read coverage
  d <- density(x$tumour, n = nrow(x))
  maxpeak <- d$x[which.max(d$y)]
  
  ## Record the bins used to compute the maxpeak to use later
  maxpeak_ind <- which(x$tumour < (maxpeak + d$bw) & x$tumour > (maxpeak - d$bw))

  ## Remove tibble formatting
  x <- as.data.frame(x)
  

  ## Get median normal coverage
  median_normal <- median(x$normal)
  
  ## Compute normalization factors by comparing the individual
  norm_factor <- x$normal/median_normal
  
  ## Correct tumour depth by the norm factor
  x$tumour <- x$tumour/norm_factor
  
  # Extract X chromosome regions
  ## TODO: Make this compatible with running without X-chromosomes
  ## Filter X chromosome data
  chrX <- x[x$chr == "X",]
  ## Decide whether this is a male or not
  isMale=F
  if(which.min(abs((median(chrX$window_size) - median(x$window_size)) - c(0, median(x$window_size)))) == 2){
    x$normal[x$chr == "X"] <- x$normal[x$chr == "X"]/2
    isMale=T
  }
  ## Print sex result
  if(verbose){
    if(isMale){
      print("Automated sex detection infers MALE sex")
    }else{
      print("Automated sex detection infers FEMALE sex")
    }
  }

  ## Rename columns to tumour/normal
  x <- x %>% dplyr::rename("y_raw" = "tumour", "x_raw" = "normal")
  
  ## Debugging plots, hardly used anymore anyways
  ## Shows the data normalization progress
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
  
  ## Corrected read depth values
  x$residual <- GCnorm$residuals
  
  ## Detect offset of residual value from previous, in case data got shifted
  d <- density(x$residual[maxpeak_ind], n = 2^16)
  offset <- d$x[which.max(d$y)]
  
  ## Holdover from old code where normalization was applied to bin sizes as well
  x$normalized_size <- x$window_size
  
  
  
  ## Another debug plot
  if(debugPlots){
    GCnorm2plot <- ggplot(x, aes(y = residual, x = normalized_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Linear normalization of read counts by GC content")
    print(GCnorm2plot)
  }
  ## Apply corrections to depth, recompute maxpeak
  x <- x %>% mutate(residual = residual + maxpeak - offset) %>% dplyr::rename("corrected_depth" = "residual")
  d <- density(x$corrected_depth, n = 2^16)
  maxpeak <- d$x[which.max(d$y)]
  output <- list("x" = x, "maxpeak" = maxpeak, "merged" = closest)
  return(output)
}
