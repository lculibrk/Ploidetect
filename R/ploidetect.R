#' @import data.table
#
#'@export
ploidetect_presegment <- function(all_data, centromeres = F, simplify_size = 100000){
  
  ## Load centromeres
  if(centromeres != F){
    centromeres = fread(centromeres)
    names(centromeres) = c("chr", "pos", "end", "band", "type")
    centromeres = centromeres[type == "acen"]
  }

  ## Add columns to input data
  names(all_data) <- c("chr", "pos", "end", "tumour", "normal", "maf", "gc")
  
  ## Run ploidetect_preprocess
  
  output <- ploidetect_preprocess(all_data = all_data, centromeres = centromeres, simplify_size = simplify_size, debugPlots = T)
  
  ## Unpack output
  processed_data <- output$x
  maxpeak <- output$maxpeak
  n_merged <- output$merged
  
  
  ## Initial prior, low TC so we get conservative segmentation
  tp=0.15
  ploidy=2
  
  ## Segment the genome!
  segmented_data <- ploidetect_segmentator(processed_data, maxpeak = maxpeak, verbose = F, tp = tp, ploidy = ploidy, segmentation_threshold = 0.25)
  segmented_data <- do.call(rbind.data.frame, segmented_data) %>% ungroup()
  return(list("segmented_data" = segmented_data, "maxpeak" = maxpeak, "all_data" = all_data))
}

#'@export
ploidetect <- function(in_list, call_cna = F){
  ## Unpack in_list
  segmented_data <- in_list$segmented_data
  
  den <- density(segmented_data$corrected_depth, n = 2^16, bw = "nrd0")
  all_data <- in_list$all_data
  
  maxpeak <- den$x[which.max(den$y)]
  
  #segmented_data %>% dplyr::filter((end - pos) < 1.5e+06, corrected_depth < 1e+05) %>% ggplot(aes(x = end - pos, y = corrected_depth)) + geom_point(size = 1, alpha = 0.1)
  
  
    
  ## Initialize plots object
  plots <- list()
  
  segments <- segmented_data %>% group_by(chr, segment) %>% dplyr::summarise("pos" = first(pos), "end" = last(end), "var_mafs" = sd(abs(as.numeric(unlist(unmerge_mafs(maf))) - 0.5) + 0.5, na.rm = T), "maf" = merge_mafs(maf, na.rm = T, exp = T), var_segments = sd(corrected_depth), segment_depth = median(corrected_depth))
  #segment_dp_sd <- segments$var_segments^2 %>% mean(na.rm=T) %>% sqrt()
  segment_maf_sd <- segments$var_mafs^2 %>% mean(na.rm = T) %>% sqrt()
  
  ## Compute limits on TP to sweep over
  ## Compute 5 and 95th percentiles
  quant <- quantile(segmented_data$corrected_depth, probs = c(0.05, 0.95))
  
  ## Compute initial constraints on tp
  # Assume pure, low ploidy tumour
  max_tp = get_coverage_characteristics(tp = 1, ploidy = 1, maxpeak = maxpeak) %>% tp_diffploidy(2)
  
  
  # Assume impure, high ploidy tumour
  # Risk here is to waste time by computing too low TP models which are unable to explain the data
  #sweep_tp = 0
  #not_identified=T
  #while(not_identified){
  #  sweep_tp = sweep_tp + 0.1
  #  min_tp_char = get_coverage_characteristics(tp = sweep_tp, ploidy = 5, maxpeak = maxpeak)
  #  min_tp_limits = min_tp_char$cn_by_depth[c(1, 11)]
  #  not_identified <- any(c(min_tp_limits[2] < quant[2], min_tp_limits[1] > quant[1]))
  #}
  
  ## Where to start the sweep
  #min_tp <- min_tp_char %>% tp_diffploidy(2)
  
  min_dp <- seq(from = quant[1], to = quant[2], length.out = 11) %>% diff() %>% mean()
  
  max_dp <- get_coverage_characteristics(tp = 1, ploidy = 1, maxpeak = maxpeak)$diff
  
  #max_dp <- quant[2] - quant[1]
  
  model_params <- data.frame("modeling_tp" = c(), "depth_lik" = c(), "maf_lik" = c(), "jitter" = c(), "tp" = c(), "ploidy" = c(), lowest = c(), "n_imputed" = c())
  pdfs <- list()
  models <- list()
  exp_plots <- list()
  metric <- c()
  mads <- c()
  n50s <- NULL
  collected_sds <- c()
  CNA_p <- list()
  
  
  if(exists("new_CN_calls")){
    rm(new_CN_calls)
  }
  segs <- unlist(lapply(split(1:nrow(segmented_data), segmented_data$chr), function(x)1:length(x)))
  segmented_data$segment_depth <- segmented_data$corrected_depth
  

  for(d_diff in seq(from = min_dp, to = max_dp, length.out = 41)){
    print(d_diff)
#  for(d_diff in out_models$modeling_tp[1:10]){
    #Removed to change sweep from "tp" to depth-based{
    ## Initialize parameters
    # Coverage characteristics
    #cov_char <- get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)
    # Predicted depths
    #predictedpositions <- cov_char$cn_by_depth
    predictedpositions <- seq(from = (maxpeak - 5 * d_diff), to = (maxpeak + 5 * d_diff), by = d_diff)
    names(predictedpositions) <- 1:11
    # Data.frame
    data_den <- density(segmented_data$corrected_depth, n = nrow(segmented_data))
    # Compute SD by comparison with KDE
    
    em_sd <- match_kde_height(data = segmented_data$corrected_depth, means = predictedpositions, sd = data_den$bw)
    
    
    resp_mat <- compute_responsibilities(segmented_data$corrected_depth, predictedpositions, em_sd)
    
    #plot_density_gmm(segmented_data$segment_depth, means = predictedpositions, weights = colSums(resp_mat, na.rm = T), sd = em_sd)
    
    
    d_diff <- gmm_em_fixed_var(data = segmented_data$corrected_depth, means = predictedpositions, var = em_sd)
    
    predictedpositions <- seq(from = (maxpeak - 5 * d_diff), to = (maxpeak + 5 * d_diff), by = d_diff)
    names(predictedpositions) <- 1:11
    
    resp_mat <- compute_responsibilities(segmented_data$corrected_depth, predictedpositions, em_sd)
    
    em_sd <- match_kde_height(data = segmented_data$corrected_depth, means = predictedpositions, sd = em_sd)
      
    ## Detect redundant models
    if(nrow(model_params) > 0){
      closest <- which.min(abs(d_diff - model_params$modeling_tp))
      if((d_diff - model_params$modeling_tp[closest])/(d_diff + model_params$modeling_tp[closest]) < 0.01){
        next
      }
    }
    # Compute likelihood values
      
    depth_posterior <- parametric_gmm_fit(segmented_data$corrected_depth, means = predictedpositions, variances = em_sd)
    

    # Compute chi-squared probabilities from chi values
    #tot_probs <- depth_posterior
    # Detect ties, where a point fits nothing
    ties <- depth_posterior %>% apply(., 1, function(x){all(x == 0)}) %>% which()
    # Create new df from original data
    plot_data <- segmented_data %>% mutate("fit_peak" = apply(depth_posterior, 1, which.max))
    # Enter ties data
    plot_data$fit_peak[ties] <- NA
    # Compute weights of the mixture model
    #proportions <- table_vec(plot_data$fit_peak)
    #props_names <- names(proportions)
    #proportions <- proportions/sum(proportions)
    
    # Compute posterior probabilities
    #tot_probs <- 1 - tot_probs
    #tot_probs <- lapply(1:nrow(depth_posterior), function(x){
    #  x <- tot_probs[x,as.numeric(props_names)]
    #  x <- ((x * proportions)/sum(x * proportions)) * x
    #}) %>% do.call(rbind, .)
    
    #tot_probs
    
    plot_data$prob <- apply(depth_posterior, 1, max)
    
    plot_data$prob[ties] <- 0
    
    plot_data$fit_peak[!ties] <- depth_posterior[!ties,] %>% as.matrix() %>% apply(., 1, which.max)
    
    probs <- plot_data$prob
    
    #residuals <- outer(plot_data$corrected_depth, predictedpositions, FUN = "-") %>% abs
    #residuals <- residuals %>% apply(1, min)
    #residuals <- data.frame("residual" = residuals, "n" = 1:length(residuals)) %>% filter(residual < d_diff)
    #fit_unparametric_gmm <- function(observations, components = 10){
      ## Initial parameters
    #  var=segment_dp_sd
    #  means=sort(observations[sample(1:length(observations), components)])
    #  weights=1
    #  ## Generate gmm fit with initial parameters
    #  fits <- parametric_gmm_fit(data = observations, means = means, variances = var)
    #  nofit <- fits %>% apply(., 1, function(x){all(x ==1)}) %>% which
    #  fits[nofit,] <- 1
    #  comp.fits <- fits %>% apply(., 1, which.min)
    #  subclonal_proportions <- table_vec(comp.fits)
    #  subclonal_proportions <- subclonal_proportions/sum(subclonal_proportions)
    #  fits <- 1 - fits
    ##  ## Estimate cumulative probability
    #  if(ncol(fits) > 1){
    #    new_fits <- lapply(1:nrow(fits), function(x){
    #      x <- fits[x,as.numeric(subclonal_proportions_names)]
    #      x <- ((x * subclonal_proportions)/sum(x * subclonal_proportions)) * x
    #    }) %>% do.call(rbind, .)
    #  }
    #  bayes_probs <- new_fits %>% apply(., 1, max)
    #  new_fits[is.na(bayes_probs),] <- 0
    #  comp.fits <- new_fits %>% apply(1, which.max)
    #  comp.fits[is.na(bayes_probs)] <- 0
    #  bayes_probs <- new_fits %>% apply(., 1, max)
    #  sum_probs <- sum(bayes_probs)
    #  df_probs <- data.frame("observations" = observations)
    #  df_probs$comp <- comp.fits
    #  means_by_comp <- df_probs %>% group_by(comp) %>% dplyr::summarise("obs" = mean(observations))
    #  means <- means_by_comp$obs
    #  names(means) <- means_by_comp$comp
    #  var_by_comp <- df_probs %>% group_by(comp) %>% dplyr::summarise("sd" = sd(observations))
    #  var <- var_by_comp$sd
    #  names(var) <- var_by_comp$comp
    #  weights <- table_vec(comp.fits)
    #  weights <- weights/sum(weights)
      
    #  print(sum_probs)
    #  
    #  pdf <- mixpdf_function(means = means, proportions = weights, sd = var)
    #  plot(pdf(min(observations):max(observations)), type = "l")
    #}
    
    
    #plot_data$fit_peak <- as.numeric(colnames(tot_probs)[plot_data$fit_peak])
    
    peak_fits <- plot_data$fit_peak
    
    #plot(density(plot_data$fit_peak, na.rm = T))
    #
    plot_data %>% ggplot(aes(x = end - pos, y = segment_depth, color = prob)) + geom_point(size = 0.1) + scale_x_continuous(limits = c(50000, 1.5e+05)) + scale_y_continuous(limits = c(0, 1.5e+05)) + scale_color_viridis()
    
    plot_data %>% ggplot(aes(x = end - pos, y = segment_depth, color = prob)) + geom_point(size = 0.1) + scale_x_continuous(limits = c(50000, 1.5e+05)) + scale_y_continuous(limits = c(0, 1.5e+05)) + scale_color_viridis() + geom_smooth(method = "lm")
    
    #plot_data %>% ggplot(aes(x = end - pos, y = segment_depth)) + geom_point(size = 0.1) + scale_x_continuous(limits = c(50000, 1.5e+05)) + scale_y_continuous(limits = c(0, 1.5e+05))
    
    peak_dp_sd <- em_sd
    
    #depth_posterior <- parametric_gmm_fit(segmented_data$segment_depth, means = predictedpositions, variances = peak_dp_sd)
    
    
    #plot_data <- plot_data[findInterval(plot_data$segment_depth, quantile(plot_data$segment_depth, probs = c(0.005, 0.995))) == 1,]
    
    #depth_posterior <- parametric_gmm_fit(plot_data$segment_depth, means = predictedpositions, variances = peak_dp_sd)
    
    
    #peak_fits <- depth_posterior^2 %>% apply(., 1, which.min)
    #plot_data$fit_peak <- peak_fits
    #ties <- depth_posterior^2 %>% pchisq(df = 2) %>% apply(., 1, function(x){all(x > 0.9)})
    #plot_data$fit_peak[ties] <- NA
    
    #peak_dp_sd <- plot_data %>% group_by(fit_peak) %>% summarise("dev" = sd(segment_depth)) %>% filter(!is.na(fit_peak) & !is.na(dev)) %>% select(dev) %>% unlist() %>% mean()
    #depth_posterior <- parametric_gmm_fit(plot_data$segment_depth, means = predictedpositions, variances = peak_dp_sd)
    #peak_fits <- depth_posterior %>% apply(., 1, function(x)which.min(x^2))
    #plot_data$fit_peak <- peak_fits
    
    
    
    
    #plot(density(peak_fits))
    #tmat <- matrix(nrow = 4, ncol = 4)
    #plot_data$prob <- depth_posterior %>% apply(., 1, function(x)min(x^2)) %>% pchisq(df = 2)
    
    #predictedpositions
    
    
    ## Compute skew for peaks > ploidy
    
    skew <- depth_posterior %>% apply(., 1, function(x)x[which.min(x^2)])
    
    sk.signs <- sign(skew)
    
    skew <- skew^2 %>% pchisq(df = 2) * sk.signs
    
    #plot(density(skew[plot_data$fit_peak > 5], na.rm = T))
    
    
    tot_skew <- mean(abs(skew[plot_data$fit_peak != 6]), na.rm = T)
    
    
    
    depth_posterior <- depth_posterior^2
    
    split_by_peak <- split(plot_data, f = plot_data$fit_peak)
    
    #mean_probs <- lapply(split_by_peak, function(x){
    #  rms(x$prob)
    #}) %>% unlist()
    
    #fit_stat <- em_result$liks
    
    resp <- depth_posterior/rowSums(depth_posterior)
    resp[is.na(resp)[,1],] <- 0
    fit_stat <- rowSums(depth_posterior * resp) %>% mean()
    
    #plot_data %>% ggplot(aes(x = end - pos, y = segment_depth, color = prob == 0)) + geom_point(size = 0.1) + scale_x_continuous(limits = c(50000, 1.5e+05)) + scale_y_continuous(limits = c(0, 1.5e+05)) + scale_color_viridis(discrete = T)
    
    prop_unfit <- (plot_data$prob < 0.05) %>% which %>% length()
    
    
    #plot(density(peak_fits))
    #proportions <- table(as.numeric(plot_data$fit_peak) - 1)
    #prop_names <- names(proportions)
    #proportions <- as.vector(proportions)
    #props <- proportions
    #prop_unfit <- prop_unfit/sum(c(proportions,prop_unfit))
    #proportions <- proportions/sum(proportions)
    #names(proportions) <- prop_names
    #proportions
    #not_in <- which(!(0:10 %in% names(proportions))) - 1
    #zer <- rep(0, times = length(not_in))
    #names(zer) <- not_in
    #proportions <- c(proportions, zer)
    #proportions <- proportions[order(as.numeric(names(proportions)))]
    #num_props <- proportions
    #proportions <- proportions/sum(proportions)
    
    #depth_posterior <- depth_posterior %>% pchisq(df = 2)
    #depth_posterior <- 1 - depth_posterior
    
    #depth_posterior <- lapply(1:nrow(depth_posterior), function(x){
    #  x <- softmax(depth_posterior[x,])
    #x <- 1 - x
    #  x <- x * proportions/sum(x * proportions)
    #}) %>% do.call(rbind, .)
    #plot_data$prob <- depth_posterior %>% apply(1, max)
    #plot_data$prob[is.na(plot_data$prob)] <- 0
    #plot_data %>% filter(fit_peak < 12) %>% ggplot(aes(x = end - pos, y = segment_depth, color = fit_peak)) + geom_point(size = 0.1) + scale_x_continuous(limits = c(50000, 1.5e+05)) + scale_y_continuous(limits = c(0, 1.5e+05)) + scale_color_viridis(discrete = F)
    
    #fit_stat <- mean(plot_data$prob)
    #num_props <- num_props[num_props != 0]
    
    ## Identify which gaussians are represented in at least 1% of DNA present
    proportions <- compute_responsibilities(segmented_data$corrected_depth, predictedpositions, em_sd)
    proportions <- colSums(proportions)/sum(colSums(proportions))
    which.1pct <- proportions > 0.01
    which.1pct <- names(proportions)[which.1pct]
    which.1pct.min <- as.numeric(which.1pct) %>% min
    which.1pct.max <- as.numeric(which.1pct) %>% max
    onepct_pos <- predictedpositions[names(predictedpositions) %in% which.1pct.min:which.1pct.max]
    onepct_weights <- proportions[names(proportions) %in% which.1pct.min:which.1pct.max]
    
    imputed_peaks <- names(onepct_pos)[!(names(onepct_pos) %in% which.1pct)]
    
    if(length(imputed_peaks) > 0){
      
      imputed_weights <- lapply(imputed_peaks, function(x){
        x <- as.numeric(x)
        mean(c(onepct_weights[which(names(onepct_weights) %in% as.character(x - 1))], onepct_weights[which(names(onepct_weights) %in% as.character(x + 1))]))
      }) %>% unlist()
      
      names(imputed_weights) <- imputed_peaks
      
      #imputed_weights <- onepct_weights[which(names(onepct_weights) %in% imputed_peaks)]
      
      n_imputed <- sum(imputed_weights)
    }else{n_imputed = 0}
    
    ploidy_fit_data <- plot_data %>% 
      dplyr::mutate(fit_peak = fit_peak, CN = fit_peak) %>% 
      dplyr::select(chr, segment, maf, fit_peak, prob) %>% 
      dplyr::filter(fit_peak %in% names(onepct_pos), !is.na(maf)) %>% 
      tidyr::separate_rows(maf, sep = ";") %>% 
      mutate(median_maf = as.numeric(maf))
    
    
    
    maf_scores <- data.frame("maf" = c(), "tp" = c(), "ploidy" = c(), "lowest" = c())
    jp_tbls <- list()
    
    
    for(lowest in 0:2){
      result <- maf_gmm_fit(depth_data = segmented_data$corrected_depth, vaf_data = segmented_data$maf, chr_vec = segmented_data$chr, means = predictedpositions, variances = em_sd, maf_variances = segment_maf_sd, lowest = lowest, maxpeak = maxpeak)
      maf_scores <- rbind.data.frame(maf_scores, result$model)
      jp_tbls[[paste0(lowest)]] <- result$jp_tbl
      next
    }
    

    ## Select best model
    maf_scores <- maf_scores[maf_scores$ploidy != 0,]
    if(nrow(maf_scores) == 0){
      next
    }
    maf_scores <- maf_scores[which.max(maf_scores$maf),]
    tp = maf_scores$tp[1]
    ploidy = maf_scores$ploidy[1]
    lowest = maf_scores$lowest[1]
    joint_probs = jp_tbls[[paste0(lowest)]]
  
    ## Generate absolute-CN positions
    #abspos <- predictedpositions[(6-ploidy):11]
    #names(abspos) <- seq(from = 0, length.out = length(abspos))
    
    ## Compute residuals of fit
    #residuals <- apply(abs(outer(segmented_data$segment_depth, abspos, "-")), 1, min)
    #out_of_range <- residuals > d_diff/2
    
    #ggplot(segmented_data, aes(x = end - pos, y = segment_depth, color = out_of_range)) + geom_point()
    
    #in_range <- 1-length(which(out_of_range))/length(residuals)
  
    #put_subclones <- residuals[!out_of_range]
    
    
    ## Traditional mixture EM
    gmm_em_residuals <- function(data, ncomps = 5){
      nsweep=data.frame("loglik" = c(), "comps" = c())
      for(i in 1:ncomps){
        initials <- kmeans(x = data, centers = i)
        ## Initialize means
        initial_means <- sort(initials$centers)
        initial_var <- aggregate(data, by = list(initials$cluster), FUN = "sd")
        initial_weights <- table_vec(initials$cluster)/sum(table_vec(initials$cluster))
        initial_var <- weighted.mean(initial_var$x, initial_weights)
        #Set variables
        vars <- initial_var
        weights <- initial_weights
        means <- initial_means
        # Perform the first iteration to get log-likelihood score
        fit <- parametric_gmm_fit(data = data, means = means, variances = vars)
        responsibilities <- fit/rowSums(fit)
        weights <- colSums(responsibilities)
        weights <- weights/sum(weights)
        log_lik <- sum(log(rowSums(t(t(fit) * weights))))
        results <- list(list("means" = means, "weights" = weights, "vars" = vars, "loglik" = log_lik))
        iters <- 1
        while(iters < 10){
          means = colSums(responsibilities * data)/colSums(responsibilities)
          vars = sqrt(colSums(responsibilities * (data - means)^2)/colSums(responsibilities))
          fit = parametric_gmm_fit(data = data, means = means, variances = vars)
          responsibilities <- fit/rowSums(fit)
          weights <- colSums(responsibilities)
          weights <- weights/sum(weights)
          log_lik <- sum(log(rowSums(t(t(fit) * weights))))
          iters <- iters+1
          results[iters] <- list(list("means" = means, "weights" = weights, "vars" = vars, "loglik" = log_lik))
        }
        f <- mixpdf_function(means = means, proportions = weights, sd = vars)
        plot(f(0:4000), type = "l"); rug(x = data, col = rgb(0,0,0,0.01))
        results <- results[[lapply(results, function(x)x$loglik) %>% unlist %>% which.max]]
        nsweep=rbind.data.frame(nsweep, data.frame("loglik" = results$loglik, "comps" = i))
      }
      nsweep$bic <- -2*nsweep$loglik + (nsweep$comps)*log(length(data))
      plot(bic ~ comps, data = nsweep)
      
      
    }
    
    ## Fix proportions
    
    toadd_names <- c(1:11)[!(1:11 %in% names(proportions))]
    toadd <- rep(0, times = length(toadd_names))
    names(toadd) <- toadd_names
    
    proportions <- c(toadd, proportions)
    proportions <- proportions[order(as.numeric(names(proportions)))]
    
    pdf_fun <- mixpdf_function(means = predictedpositions, proportions = proportions, sd = peak_dp_sd)
    #cdf_fun <- function(x){
    #  return(cumsum(pdf_fun(x)))
    #}
    data_den <- density(segmented_data$corrected_depth, n = nrow(segmented_data))
    #plot(data_den)
    data_den <- data.frame("x" = data_den$x, "data" = data_den$y)
    #data_den$pdf <- pdf_fun(data_den$x)
    #fit_stat <- cosine_sim(data_den$data, data_den$pdf)
    #fit_stat <- bhatt_dist(data_den$data, data_den$pdf)
    #fit_stat <- depth_posterior %>% apply(1, min) %>% pchisq(df = 2) %>% rms()
    #data_den$diff <- abs(data_den$data - data_den$pdf)
    #data_den$diff <- apply(data_den, 1, function(x){x[4]/max(x[2], x[3])})
    
    #positions_indices <- sapply(predictedpositions, function(x) data_den$diff[which.min(abs(data_den$x - x))])
    
    #pdf_fun <- mixpdf_function(means = predictedpositions, proportions = 1, sd = peak_dp_sd)
    data_den$pdf <- pdf_fun(data_den$x)$y
    

    
    plot_posmafs <- plot_data %>% group_by(fit_peak) %>% dplyr::summarise(maf = round(mean(as.numeric(unlist(unmerge_mafs(maf, flip = T)))), 2), med_cov = median(segment_depth)) %>% mutate("pdf_val" = pdf_fun(med_cov)$y) %>% filter(!is.na(maf))
    plot_posmafs$dat <- 0
    for(i in 1:nrow(plot_posmafs)){
      plot_posmafs$dat[i] <- data_den$data[which.min(abs(data_den$x - plot_posmafs$med_cov[i]))]
    }
    
    dodge_amount <- quantile(data_den$data, c(0, 1))[2]/25

    

    
    
    
    ## Try segmentation for evaluation
    unaltered <- round(mean(segmented_data$end - segmented_data$pos)/mean(all_data$end - all_data$pos), digits = 0)
    
    diffsum <- function(x){
      n <- c()
      for(i in 2:length(x)){
        n[i-1] <- x[i] + x[i-1]
      }
      n
    }
    diffmin <- function(x){
      n <- c()
      for(i in 2:length(x)){
        n[i-1] <- min(c(x[i], x[i-1]))
      }
      n
    }
    
    #new_CN_calls$segment_depth
    
    #new_CN_calls %>% group_by(CN, segment) %>% dplyr::summarise(dev = sd(corrected_depth), n = n()) %>% ungroup %>% dplyr::select(3,4) %>% summarise(weighted.mean(dev, w = n, na.rm = T))
    
    data_den <- density(segmented_data$corrected_depth, n = nrow(segmented_data))
    #plot(data_den)
    data_den <- data.frame("x" = data_den$x, "data" = data_den$y)
    #data_den$pdf <- pdf_fun(data_den$x)
    #fit_stat <- cosine_sim(data_den$data, data_den$pdf)
    #fit_stat <- bhatt_dist(data_den$data, data_den$pdf)
    #fit_stat <- depth_posterior %>% apply(1, min) %>% pchisq(df = 2) %>% rms()
    #data_den$diff <- abs(data_den$data - data_den$pdf)
    #data_den$diff <- apply(data_den, 1, function(x){x[4]/max(x[2], x[3])})
    
    #positions_indices <- sapply(predictedpositions, function(x) data_den$diff[which.min(abs(data_den$x - x))])
    
    #pdf_fun <- mixpdf_function(means = predictedpositions, proportions = 1, sd = peak_dp_sd)
    data_den$pdf <- pdf_fun(data_den$x)$y
    
    p <- data_den %>% 
      filter(x < max(predictedpositions) + predictedpositions[2] - predictedpositions[1]) %>% 
      ggplot() + 
      geom_line(aes(x = x, y = pdf, color = "predicted distribution")) + 
      geom_line(aes(x = x, y = data)) + 
      geom_vline(xintercept = predictedpositions, linetype = 2, alpha = 0.2) + 
      geom_text(data = plot_posmafs, aes(x = med_cov, y = dat + dodge_amount, label = paste0("MAF = ", maf))) + 
      scale_x_continuous(limits = c(min(data_den$x)-1, max(predictedpositions) + predictedpositions[2] - predictedpositions[1])) + 
      ggtitle(paste0("TP = ", round(maf_scores$tp[1], digits = 2), ", Ploidy = ", maf_scores$ploidy[1])) + 
      xlab("Read Depth") + 
      ylab("Density") + 
      theme_cowplot() + 
      theme(axis.text.x = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title = element_text(size=20),
            plot.title = element_text(size = 15),
            legend.position = "none")

    
    print(p)
    
    train_df <- data.frame("segment_depth" = predictedpositions, "CN" = seq(from = ploidy - 5, to = ploidy + 5))
    
    model <- lm(CN ~ segment_depth, train_df)
    
    seg_mat <- joint_probs
    seg_mat <- seg_mat/rowSums(seg_mat)
    seg_mat <- as.matrix(seg_mat)
    seg_mat[is.na(seg_mat)] <- 0
    seg_mat <- data.table(seg_mat)
    
    
    
    segs <- ploidetect_prob_segmentator(prob_mat = seg_mat, ploidy = ploidy, chr_vec = segmented_data$chr, seg_vec = 1:nrow(segmented_data), dist_vec = segmented_data$corrected_depth, subclones_discovered = F, lik_shift = 1.5)
    
    new_seg_data <- segmented_data %>% mutate("segment" = segs, "prob" = probs) %>% group_by(chr, segment) %>% dplyr::mutate(segment_depth = median(corrected_depth), prob = median(prob))
    
    new_seg_data %>% ggplot(aes(x = pos, y = segment_depth)) + geom_point() + facet_wrap(~chr)
    
    practical_sd <- new_seg_data %>% group_by(chr, segment) %>% add_tally() %>% filter(n > 2) %>% dplyr::summarise("var_dp" = median(abs(corrected_depth - median(corrected_depth))), "n" = first(n)) %>% ungroup %>% dplyr::summarise("var_dp" = weighted.mean(var_dp, w = n)) %>% unlist()
    mad_gmm <- gmm_mad(new_seg_data$corrected_depth, means = predictedpositions, variances = em_sd)
    
    new_seg_data %>% filter(chr == "3") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis()
    
    contiguity_dat <- data.table(new_seg_data)[, .(pos = first(pos), end = last(end)), by = list(chr, segment)]
    
    n50 <- nx_fun(as.numeric(contiguity_dat$end - contiguity_dat$pos), 0.5)
    
    if(exists("n50s")){
      n50s <- c(n50s, n50)
    }else{
      n50s <- n50
    }
    #new_seg_data %>% group_by(chr, segment) %>% add_tally() %>% filter(n > 2) %>% dplyr::summarise("var" = sd(corrected_depth), n = first(n))
    
    if(exists("mads")){
      mads <- c(mads, mad_gmm)
    }else{
      mads <- mad_gmm
    }
    if(exists("collected_sds")){
      collected_sds <- c(collected_sds, practical_sd)
    }else{
      collected_sds <- practical_sd
    }
    
    exp_plots <- c(exp_plots, list(p))
    
    #diff_df <- data.frame("CN_d" = diff(seg_diff$CN_r), "dp" = diff(seg_diff$dp), "n" = diffmin(seg_diff$n)) %>% mutate("cor_dp" = dp/CN_d)
    
    #cn_diff <- weighted.mean(diff_df$cor_dp, w = diff_df$n)
    
    #agreement <- abs(cn_diff - d_diff)/(d_diff)
    
    
    #predictedpositions
    #compute_normals_overlap(x = data_den$x, mu1 = predictedpositions[6], mu2 = predictedpositions[7], w1 = proportions[6], w2 = proportions[7], sd1 = peak_dp_sd, sd2 = peak_dp_sd)
    #mat <- pairwise_overlaps(x = data_den$x, positions = predictedpositions, weights = proportions, sd = peak_dp_sd)
    #overlaps <- apply(mat, 1, max, na.rm = T)
    #pheatmap(mat, cluster_cols = F, cluster_rows = F)
    ##pdf_fun(1:200000) %>% data.frame("x" = 1:length(.), "probs" = .) %>% ggplot(aes(x = x, y = probs)) + geom_line(aes(color = "predicted distribution")) + geom_vline(xintercept = predictedpositions)
    ##fit_stat <- ks.test(segmented_data$segment_depth, cdf_fun)$statistic
    ##plot_data %>% ggplot(aes(x = end - pos, y = segment_depth, color = factor(fit_peak - 1))) + geom_point(size = 0.1) + scale_x_continuous(limits = c(0, 150000)) + scale_color_viridis(discrete = T) + scale_y_continuous(limits = c(50000, 200000))
    ##depth_posterior <- depth_posterior[-which(ties),]
    ##depth_posterior <- depth_posterior[-ties,]
    ##segmented_data$segment_depth[ties]
    ##depth_probs <- depth_posterior %>% apply(., MARGIN = 1, min) %>% pchisq(df = 11) %>% mean()
    ##cn_assignments = max.col(-depth_posterior) - 1
    ##fit_stat <- sum(log(props)) - 2*log(fit_stat)
    
    comp_resps <- colSums(compute_responsibilities(new_seg_data$segment_depth, predictedpositions, em_sd))
    comp_resps <- round(comp_resps/sum(comp_resps), 5)
    #abs(sort(desc(comp_resps)))
    if(abs(sort(desc(comp_resps)))[2] < 0.01){
      model_params <- rbind.data.frame(model_params, data.frame("modeling_tp" = NA, "depth_lik" = fit_stat, "maf_lik" = maf_scores$maf[1], "jitter" = em_sd, "tp" = maf_scores$tp[1], "ploidy" = maf_scores$ploidy[1], "n_imputed" = n_imputed))
      next
    }
    
    model_params <- rbind.data.frame(model_params, data.frame("modeling_tp" = d_diff, "depth_lik" = fit_stat, "maf_lik" = maf_scores$maf[1], "jitter" = em_sd, "tp" = maf_scores$tp[1], "ploidy" = maf_scores$ploidy[1], "n_imputed" = n_imputed))
  }
  
  #candidate_models <- model_params %>% group_by(npeaks) %>% summarise(rank = which.min(depth_lik), tp = tp[rank], depth_lik = depth_lik[rank], prop_unfit = prop_unfit[rank], jitter = jitter[rank], mean_bleed = mean_bleed[rank]) %>% select(-rank)
  metric <- abs(collected_sds - mads)
  metric <- metric/max(metric)
  
  
  max_n50 <- data.table(new_seg_data)[, .(pos = first(pos), end = last(end)), by = list(chr, segment)]
  
  max_n50 <- nx_fun(sort(max_n50$end - max_n50$pos), n = 0.5)
  
  max_n50 <- max(n50s)

  
  out_models <- model_params %>%  mutate(model_stat = (depth_lik * maf_lik * (1-(mads/mad(new_seg_data$corrected_depth))) * n50s/max_n50) * 1/(1 + n_imputed), order = 1:nrow(.)) %>% bind_cols(data.frame("seg_agreement" = metric)) %>% filter(tp > 0.05) %>% arrange(desc(model_stat)) %>% mutate(model_stat = model_stat/sum(model_stat)) 
  
  out_models$model_stat <- out_models$model_stat/sum(out_models$model_stat)
  
  exp_plots <- exp_plots[out_models$order]
  
  
  
  if(is.na(out_models$modeling_tp[1])){
    out_models <- "No CNVs detected. Unable to estimate copy number"
  }
  

  #cna_calls <- ploidetect_cna_sc(all_data = in_list$all_data, segmented_data = segmented_data, tp = out_models$tp[1], ploidy = out_models$ploidy[1], maxpeak = maxpeak)
  
  return(list("model_info" = out_models, "plots" = exp_plots, "maxpeak" = maxpeak, "segmented_data" = segmented_data))
}
