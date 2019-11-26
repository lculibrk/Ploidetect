#' Z-score
#' 
#' \code{zt} returns the z-score (number of standard deviations from the mean)
#' of the input given a mean and variance
#' @param obs (numeric) value or vector of values
#' @param mean (numeric) the mean that should be used for calculation of z
#' @param sd (numeric) standard deviation for calculation of z
#' @param absval (boolean) should the absolute value of the z-score be returned?
#' @return (numeric) a vector of equal length to the input
#' @examples 
#' zt(1, 1, 1)
#' zt(0, 1, 1)
#' zt(10, 1, 1)
#' zt(1:5, 1, 0.3)
zt <- function(obs, mean, sd, absval = T){
  if(absval){
    abs((obs-mean)/sd)
  }else{
    (obs-mean)/sd
  }
}
#' Standardized normal p-value
#' 
#' \code{zt_p} returns \eqn{p(x)} by mapping x to a normal distribution with mean
#' equal to zero and standard deviation equal to one using its z-score
#' @param obs (numeric) value or vector of values
#' @param mean (numeric) the mean of the gaussian being compared against
#' @param sd (numeric) standard deviation of the gaussian being compared
#' against
#' @return (numeric) a p-value between zero and ~0.399
#' @examples 
#' zt_p(1, 1, 1)
#' zt_p(0, 1, 1)
#' zt_p(1:5, 1, 1)
zt_p <- function(obs, mean, sd){
  dnorm(zt(obs, mean, sd), mean = 0, sd = 1)
}

#' Likelihood function of a gaussian distribution
#' 
#' \code{norm_lik} returns the value of \eqn{L(x | \mu, \sigma)}, the
#' likelihood function of a value given a gaussian distribution
#' @param all_obs (numeric) vector of values over which to compute the 
#' likelihood
#' @param mean (numeric) mean of the gaussian
#' @param sd (numeric) standard deviation of the gaussian
#' @return (numeric) vector of values equal in length to x
norm_lik <- function(all_obs, mean, sd){
  ((2*pi*sd^2)^(-length(all_obs)/2))*exp(-1/(2*sd^2)*sum(all_obs - mean))
}

#' Probability density function of a gaussian distribution
#' 
#' \code{normpdf} returns the value of \eqn{f(x, \mu, \sigma^2)}, the
#' probability density function (PDF) of a gaussian distribution
#' @param x (numeric) vector of values over which to compute the PDF
#' @param mean (numeric) mean of the gaussian
#' @param sd (numeric) standard deviation of the gaussian
#' @return (numeric) vector of values equal in length to x
normpdf <- function(x, mean, sd){
  prob <- 1/sqrt(2 * pi * sd^2) * exp(-((x - mean)^2)/(2*sd^2))
  return(prob)
}

#' Fits a gaussian mixture model with specified parameters to the data agnostic
#' to weights
#' \code{parametric_gmm_fit} returns a matrix of standardized normal 
#' probabilities with rows as observations and columns as mixture components
#' @param data (numeric) vector of values to compute probailities for
#' @param means (numeric) vector of gaussian means, one per component
#' @param variances (numeric) vector of gaussian variances, either length one
#' or equal to the length of means
#' @return matrix of values with row count equal to the length of data and
#' column count equal to the length of the means
#' @examples 
#' parametric_gmm_fit(0:5, c(1, 5), 2)
parametric_gmm_fit <- function(data, means, variances){
  gaussians <- cbind(means, variances)
  result <- outer(data, gaussians[,1], "zt_p", "sd" = gaussians[,2])
  return(result)
}
#' Computes weighted median absolute deviations (MADs) for each component of a
#' gaussian mixture model
#' 
#' \code{gmm_mad} returns the weighted median absolute deviation of data fit to
#' the specified mixture model
#' @param data (numeric) vector of values to compute the MAD over
#' @param means (numeric) vector of means for the mixture moel
#' @param variances (numeric) vector of gaussian variances, either length one
#' or equal to the length of means
#' @return a single numeric of the MAD for the model
#' @examples 
#' gmm_mad(0:5, c(1, 5), 2)
gmm_mad <- function(data, means, variances){
  responsibilities <- compute_responsibilities(data, means, variances)
  weights <- colSums(responsibilities)/sum(colSums(responsibilities))
  devs <- abs(outer(data, means, "-"))
  mad <- lapply(1:ncol(devs), function(x)devs[,x]*responsibilities[,x]) %>% do.call(cbind, .) %>% rowSums %>% median
  return(mad)
}
#' Optimizes mean differential for a special case of gaussian mixture models
#' 
#' (GMMs) where component means differ by a constant factor
#' \code{gmm_em_fixed_var} returns a difference optimized using EM and a grid
#' search
#' @param data (numeric) vector of values to for the GMM to be fit to
#' @param means (numeric) vector of means for GMM components
#' @param var (numeric) vector of gaussian variances, either length one
#' or equal to the length of means
#' @return a single numeric value representing the difference in means
gmm_em_fixed_var <- function(data, means, var){
  i_diff <- mean(diff(means))
  p_mat <- parametric_gmm_fit(data, means, var)
  r_mat <- p_mat/rowSums(p_mat)
  
  lik <- mean(rowSums(p_mat * r_mat, na.rm = T))
  
  n_means = colSums(r_mat * data, na.rm = T)/colSums(r_mat, na.rm = T)
  n_diff <- mean(diff(n_means), na.rm = T)
  n_means = seq(from = means[6] - n_diff*5, to = means[6] + n_diff*5, by = n_diff)
  p_mat <- parametric_gmm_fit(data, n_means, var)
  r_mat <- p_mat/rowSums(p_mat)
  
  n_lik <- mean(rowSums(r_mat * p_mat, na.rm = T))
  
  while(n_lik > lik){
    lik <- n_lik
    n_means = colSums(r_mat * data, na.rm = T)/colSums(r_mat, na.rm = T)
    n_diff <- mean(diff(n_means), na.rm = T)
    n_means = seq(from = means[6] - n_diff*5, to = means[6] + n_diff*5, by = n_diff)
    p_mat <- parametric_gmm_fit(data, n_means, var)
    r_mat <- p_mat/rowSums(p_mat, na.rm = T)
    n_lik <- mean(rowSums(r_mat * p_mat, na.rm = T))
  }
  condition = T
  t <- c()
  for(i in seq(from = i_diff, to = n_diff, length.out = 20)){
    p_mat <- parametric_gmm_fit(data, seq(from = means[6] - i*5, to = means[6] + i*5, by = i), var)
    r_mat <- p_mat/rowSums(p_mat)
    lik = mean(rowSums(p_mat * r_mat, na.rm = T))
    t <- c(t, lik)
  }
  out_diff <- seq(from = i_diff, to = n_diff, length.out = 20)[which.max(t)]
  return(out_diff)
}
 #' Computes responsibilities for a gaussian mixture model with specified
 #' parameters
 #' 
 #' \code{compute_responsibilities} returns a matrix as in parametric_gmm_fit
 #' except values are responsibilities rather than probabilities
 #' @inheritParams parametric_gmm_fit
 #' @return matrix of values with row count equal to the length of data and
 #' column count equal to the length of the means
compute_responsibilities <- function(data, means, variances){
  mat <- parametric_gmm_fit(data = data, means = means, variances = variances)
  resp <- mat/rowSums(mat)
  resp[is.na(resp)[,1],] <- 0
  return(resp)
}

#' Determines variance for a gaussian mixture model (GMM) given fixed means 
#' using a variation of expectation-maximization for the special case of all
#' equal variances for all components
#' 
#' \code{gmm_em_only_var} returns a value of standard deviation for the given
#' GMM
#' @param data (numeric) vector of values to for the GMM to be fit to
#' @param means (numeric) vector of means for GMM components
#' @param var (numeric) single value, initial estimate for variance
#' @return a single numeric value representing the difference in means
gmm_em_only_var <- function(data, means, initial_var){
  ## Record initial value for variances
  initial_vars <- initial_var
  initial_weights <- rep(1, times = length(means))
  #
  means <- means
  weights <- initial_weights
  vars <- initial_vars
  #
  old_lik <- 0
  new_lik <- 0.002
  liks <- data.frame("vars" = c(), "lik" = c())
  #
  iter <- 0
  while((((new_lik - old_lik)/old_lik)) > 0.001 | iter < 2){
    old_lik <- new_lik
    ## LIKELIHOOD
    gaussian_probs <- parametric_gmm_fit(data = data, means = means, variances = vars)

    ## EXPECTATION STEP
    ## Compute responsibilities

    responsibilities <- t(t(gaussian_probs) * weights)
    responsibilities <- responsibilities/rowSums(responsibilities)
    responsibilities[is.na(responsibilities)[,1],] <- 0

    new_lik <- rowMeans(gaussian_probs * responsibilities) %>% mean()
    ### Record parameter information
    liks <- rbind.data.frame(liks, data.frame("vars" = vars, "liks" = new_lik))
    
    if(all(gaussian_probs == 0)){
      vars = vars * 2
      next
    }
    
    ## MAXIMIZATION STEP
    ## Compute mixing proportions
    resp_sums <- colSums(responsibilities)
    resp_means <- colMeans(responsibilities)
    weights <- resp_sums/sum(resp_sums)
    ## Compute means, disabled atm
    #resp_means <- colSums(responsibilities * mafs)
    #means <- resp_means/colSums(responsibilities)
    ## Compute vars
    resp_vars <- outer(X = data, Y = means, FUN = "-")^2 * responsibilities
    resp_vars <- colSums(resp_vars)
    resp_vars <- sqrt(resp_vars/(resp_sums - resp_means))
    vars <- resp_vars * weights
    vars <- weighted.mean(vars, w = weights)
    iter <- iter + 1
  }
  pdf_fun <- mixpdf_function(means = means, proportions = weights, sd = vars)
  return(liks[nrow(liks) - 1,])
  plot(pdf_fun(seq(from = 0, to = 1, by = 0.01)), type = "l")
  rug(x = mafs)
}
#' Generate the probability density function (pdf) of the specified gaussian
#' mixture model (GMM)
#' 
#' \code{mixpdf} returns a vector of probability density over the values of x
#' @param x (numeric) vector of values over which to compute the pdf
#' @param means (numeric) vector of means for GMM components
#' @param proportions (numeric) vector of GMM component weights
#' @param sd (numeric) vector or single value of GMM component variance
#' @return a single numeric value representing the difference in means
mixpdf <- function(x, means, proportions, sd){
  parameters <- cbind.data.frame(means, proportions)
  parameters <- split(parameters, 1:nrow(parameters))
  out <- lapply(parameters, function(params){
    normpdf(x, params$means, sd)*params$proportions
  })
  out <- rowSums(rbindlist(cbind, out))
  return(out)
}
#' Generate a function that outputs the probability density function of a GMM
#' given an input range
#' 
#' \code{mixpdf_function} returns a vector of probability density over the values of x
#' @param means (numeric) vector of means for GMM components
#' @param proportions (numeric) vector of GMM component weights
#' @param sd (numeric) vector or single value of GMM component variance
#' @return a single numeric value representing the difference in means
mixpdf_function <- function(means, proportions, sd){
  outfun <- function(x){
    parameters <- cbind.data.frame(means, proportions, sd)
    parameters <- split(parameters, 1:nrow(parameters))
    out <- lapply(parameters, function(params){
      normpdf(x, params$means, params$sd)*params$proportions
    })
    out <- do.call(cbind, out) %>% apply(1, sum)
    out=data.frame("x" = x, "y" = out)
    return(out)
  }
  return(outfun)
}
#' Compute the root-mean-squared (rms) deviation of a vector
#' 
#' \code{rms} returns the rms of the input vector
#' @param x (numeric) input vector
#' @param na.rm (boolean) should NAs be removed?
#' @return the rmsd of the input
rms <- function(x, na.rm = T){
  sqrt(mean(x^2, na.rm = na.rm))
}
#' Convert vector into probability distribution
#' 
#' \code{props} returns the proportion-transformed version of input vector
#' @param x (numeric) input vector
#' @return a numeric vector of equal length to input
props <- function(x){
  x/sum(x)
}

#' Compute predicted depths based on coverage characteristics
#' 
#' \code{depth} returns the depth of of each specified copy number given the 
#' input parameters
#' @param maxpeak the depth of the most common copy number of the genome 
#' (eg. where CN = Ploidy)
#' @param d the depth difference per copy number. Can be calculated from 
#' \code{get_coverage_characteristics}
#' @param P the ploidy (most common copy number)
#' @param n a vector of all copy numbers to compute depths for
#' @return a numeric vector of predicted depths for each n, with names equal to 
#' n
depth <- function(maxpeak, d, P, n){
  return(structure(maxpeak - d * (P - n), names = n))
}

#' Compute tumour purity at a different ploidy
#' 
#' \code{tp_diffploidy} returns the tumour purity at the given ploidy given 
#' input coverage characteristics, as returned by get_coverage_characteristics
#' @param cov_char coverage chracteristics, as returned by 
#' get_coverage_characteristics
#' @param new_ploidy new ploidy to compute tumour purity for
#' @return single numeric value for tumour purity
tp_diffploidy <- function(cov_char, new_ploidy){
  d0 = depth(maxpeak = cov_char$maxpeak, d = cov_char$diff, P = new_ploidy, n = 0)
  d2 = depth(maxpeak = cov_char$maxpeak, d = cov_char$diff, P = new_ploidy, n = 2)
  tp = d0/(d2)
  return(1 - tp)
}

#' Compute cosine similarity of two vectors
#' 
#' \code{cosine_sim} returns the cosine similarity of two vectors
#' @param vec1 the first vector
#' @param vec2 the second vector
#' @return single numeric value for cosine similarity
cosine_sim <- function(vec1, vec2){
  if(length(vec1) != length(vec2)){
    stop("Two vectors have unequal length")
  }
  a <- sum(vec1 * vec2)
  b <- sqrt(sum(vec1^2))
  c <- sqrt(sum(vec2^2))
  csim <- a/(b*c)
  return(csim)
}

#' Return a vectorized counts table
#' 
#' \code{table_vec} returns a vector of counts of unique observations of input
#' @param x the input vector to count unique occurrances of
#' @return a numeric vector of counts, where names are the counted observation
#' and values are counts
table_vec <- function(x){
  tab <- table(x)
  vec <- as.vector(tab)
  names(vec) <- names(tab)
  return(vec)
}

#' Compute the n50 length of vector of lengths
#' 
#' \code{n50_segments} returns the n50 length of input lengths
#' @param lengths input vector of lengths
#' @return the n50 value
n50_segments <- function(lengths){
  lengths=sort(lengths)
  sumlengths=sum(as.numeric(lengths))
  currentsum=0
  index=1
  while(currentsum < 0.5*sumlengths){
    currentsum = currentsum + lengths[index]
    index=index+1
  }
  return(currentsum/(index-1))
}

#' Estimate variance by comparison to a density estimate
#' 
#' \code{match_kde_height} returns standard deviation such that the probability
#' density function of the gaussian mixture model (GMM) specified by the 
#' "means" parameter closely matches the kernel density estimation for the 
#' input data
#' @param data input data to estimate density and mix the mixture model to
#' @param means the means of the GMM
#' @param sd the initial estimate of standard deviation
#' @param comparison_point value to compute density function and density 
#' estimator height for optimization
#' @return an estimation of standard deviation
match_kde_height <- function(data, means, sd, comparison_point = NA){
  ## Refine SD by comparing with KDE of data
  if(length(means) > 1){
    kde_pad <- mean(diff(means))
  }else if(length(means) == 1){
    kde_pad <- quantile(data, 0.999) - means
  }
  kde_data <- density(data[data < (max(means))], n = 2^16)
  
  ##Compute weights
  resp <- compute_responsibilities(data, means = means, variances = sd)
  proportions <- colSums(resp, na.rm = T)/sum(colSums(resp, na.rm = T), na.rm = T)
  
  if(!is.na(comparison_point)){
    maxpeak = comparison_point
  }else{maxpeak = kde_data$x[which.max(kde_data$y)]}
  
  pdf_fun <- mixpdf_function(means = means, proportions = proportions, sd = sd)
  
  plot_density_gmm(data, means, proportions, sd)
  
  ## PDF function probs at maxpeak position
  
  pdf_means <- pdf_fun(maxpeak)$y
  
  ## Closest KDE points to maxpeak
  comp_height <- kde_data$y[which.min(abs(kde_data$x - maxpeak))]
  condition <- T
  sd_min <- sd/2
  sd_max = 5 * sd
  #sd_min = sd/5
  diffs <- c()
  values_vec <- seq(from = sd_min, to = sd_max, length.out = 10)
  values_diff <- diff(values_vec)[1]
  ## rearrange values_vec
  tmp_vec <- values_vec
  values_vec <- c()
  while(length(tmp_vec) > 0){
    values_vec <- c(values_vec, tmp_vec[1])
    values_vec <- c(values_vec, tmp_vec[length(tmp_vec)])
    tmp_vec <- tmp_vec[-c(1, length(tmp_vec))]
  }
  i = 1
  while(i <= length(values_vec)){
    resp <- compute_responsibilities(data, means = means, variances = values_vec[i])
    proportions <- colSums(resp)/sum(colSums(resp))
    pdf_fun <- mixpdf_function(means = means, proportions = proportions, sd = values_vec[i])
    ## PDF function probs at means
    pdf_means <- pdf_fun(maxpeak)$y
    ## Closest KDE points to maxpeak
    comp_height <- kde_data$y[which.min(abs(kde_data$x - maxpeak))]
    diffs <- c(diffs, (pdf_means - comp_height)/(pdf_means + comp_height))
    i = i+1
  }
  diffs <- diffs[order(values_vec)]
  values_vec <- values_vec[order(values_vec)]
  
  if((abs(sum(sign(diff(diffs)))) != length(diffs) - 1) & (all(diffs < 0) | all(diffs > 0))){
    return(values_vec[which.min(abs(diffs))])
  }
  if(all(diffs < 0) | all(diffs > 0)){
    return(match_kde_height(data = data, means = means, sd = values_vec[which.min(abs(diffs))]))
  }
  inflection_points <- which(diff(sign(diffs)) != 0)
  min_diff <- inflection_points[which.min(abs(diffs)[inflection_points])]
  neighbours <- c(min_diff-1,min_diff+1)
  neighbours <- neighbours[neighbours > 0]
  comp <- neighbours[which(sign(diffs[neighbours]) != sign(diffs[min_diff]))]
  if(length(comp) > 1){
    comp <- comp[which.min(abs(diffs[comp]))]
  }
  search_diffs <- sort(c(min_diff, comp))
  
  coef <- 0.5
  coef_add <- 0.25
  condition = T
  while(condition){
    new_var <- (values_vec[search_diffs[1]] * coef) + (values_vec[search_diffs[2]] * (1 - coef))
    ##
    resp <- compute_responsibilities(data, means = means, variances = new_var)
    proportions <- colSums(resp)/sum(colSums(resp))
    ##
    pdf_fun <- mixpdf_function(means = means, proportions = proportions, sd = new_var)
    
    plot_density_gmm(data, means, proportions, new_var)
    
    ## PDF function probs at means
    pdf_means <- pdf_fun(maxpeak)$y
    ## Closest KDE points to maxpeak
    comp_height <- kde_data$y[which.min(abs(kde_data$x - maxpeak))]
    
    diff <- (pdf_means - comp_height)/(pdf_means + comp_height)
    ## Decide which direction to go
    dir <- which(sign(diffs[search_diffs]) != sign(diff))
    
    if(dir == 1){
      coef <- coef + coef_add
    }else{
      coef <- coef - coef_add
    }
    coef_add <- coef_add/2
    if(abs(diff) < 0.01){
      condition <- F
    }
  }
  out_sd = new_var
  return(out_sd)
}

#' Fits a gaussian mixture-of-alleles model
#' 
#' \code{maf_gmm_fit} fits a gaussian mixture model for each possible
#' copy number state to VAF data, computes the probability of each copy number
#' and returns a matrix of probabilities for each copy number state
#' @param depth_data input read depth data for all bins
#' @param vaf_data input vaf data for all bins
#' @param chr_vec vector of chromosomes for all bins
#' @param means vector of mixture model means
#' @param variances standard deviation for the mixture model
#' @param vaf_variances standard deviation for VAF values
#' @param lowest optional unless ploidy specified; what is the copy number of 
#' the first element of means?
#' @param ploidy optional unless lowest specified; what is the most common 
#' copy number state?
#' @return a matrix of likelihoods, where columns are copy number states and 
#' rows are input observations
maf_gmm_fit <- function(depth_data, vaf_data, chr_vec, means, variances, maf_variances, maxpeak, lowest = NA, ploidy = NA){
  if(is.na(lowest) & is.na(ploidy)){
    stop("Must specify either lowest or ploidy")
  }
  # First we want to limit our VAF analysis to well-represented copy number 
  # states in the genome, so we filter out any GMM component that explains 
  # under 1% of the overall data
  ## Compute weights of initial fit from the responsibilities
  proportions = compute_responsibilities(data = depth_data, means = means, variances = variances)
  proportions = colSums(proportions)/sum(colSums(proportions))
  ## Identify components with at least 1% of weight to use as boundaries for 
  ## VAF fitting
  which.1pct <- proportions > 0.01
  which.1pct <- names(proportions)[which.1pct]
  which.1pct.min <- min(as.numeric(which.1pct))
  which.1pct.max <- max(as.numeric(which.1pct))
  onepct_pos <- means[names(means) %in% which.1pct.min:which.1pct.max]
  onepct_weights <- proportions[names(proportions) %in% which.1pct.min:which.1pct.max]
  ## Get depth difference per copy number
  d_diff <- diff(means)[1]
  # Next we identify observations that belong to the components with >1% 
  # responsibility
  ## Fit the GMM and obtain rough estimates of component membership
  peak_fits <- apply(parametric_gmm_fit(data = depth_data, means = means, variances = variances), 1, which.max)
  peak_in <- peak_fits %in% names(onepct_weights)
  ploidy_data <- data.frame("depth" = depth_data, "vaf" = vaf_data, CN = peak_fits)

  # Since the names/indices of "means" don't necessarily mean anything we need to 
  # shift the values for peak_fits to align with their assumed copy number 
  ploidy_data$CN <- factor(ploidy_data$CN)
  lowest_pk <- as.numeric(names(onepct_pos)[1])
  scale_fct <- lowest - lowest_pk
  ## If ploidy is specified then we have to do things slightly differently
  ## We assume that names(means)[1] == 1
  if(!is.na(ploidy)){
    scale_fct <- as.numeric(names(proportions[1])) - 1
    names(means) <- as.numeric(names(means)) - scale_fct
    lowest = as.numeric(names(proportions))[1]
    onepct_pos <- means
  }
  levels(ploidy_data$CN) <- as.numeric(levels(ploidy_data$CN)) + scale_fct
  names(means) <- as.numeric(names(means)) + scale_fct
  ploidy_data$CN <- as.numeric(as.character(ploidy_data$CN))
  ## Determine tp and ploidy
  test_ploidy <- as.numeric(names(means)[which.min(abs(means - maxpeak))])
  cov_char <- list(maxpeak = maxpeak, diff = d_diff)
  test_tp <- tp_diffploidy(cov_char, new_ploidy = test_ploidy)
  ## If the model is outrageous, return something indicating this
  if(test_tp > 1.05 | test_tp < 0 | is.na(test_tp)){
    return(data.frame("maf" = NULL, "tp" = NULL, "ploidy" = NULL, "lowest" = NULL))
  }
  # Next we prepare VAF data for fitting by first extracting VAF data
  maf_ind <- 1:nrow(ploidy_data)
  maf_ind <- maf_ind[!is.na(ploidy_data$vaf)]
  ## Gneerates a list where each element is a vector of VAFs that came from
  ## one bin
  maf_list <- lapply(unmerge_mafs(as.character(ploidy_data$vaf[!is.na(ploidy_data$vaf)])), as.numeric)
  ## Set names of the list to be their row index in the input and convert to a data.frame
  names(maf_list) <- maf_ind
  maf_df <- stack(maf_list)
  maf_df$ind <- as.numeric(levels(maf_df$ind)[maf_df$ind])
  ## Generate a linear model to determine general copy number from depth
  lm_df <- data.frame("CN" = as.numeric(names(means)), "depth" = means)
  model = lm(CN ~ depth, lm_df)
  cns <- table_vec(round(predict(model, data.frame("depth" = depth_data))))
  cns <- cns[as.numeric(names(cns)) >= 0]
  # Filters CN list to the nearest "sequential" CNs in the genome. Ie. the max
  # CN we consider is the highest CN that occurs twice in a row along a 
  # chromosome
  for(cn in rev(as.numeric(names(cns)))){
    vals <- which(round(predict(model, data.frame("depth" = depth_data))) == as.numeric(cn))
    if(length(vals) > 1){
      if(min(diff(vals)) == 1){
        cns <- cns[as.numeric(names(cns)) <= cn]
        break
      }
    }
  }
  max_cns <- max(as.numeric(names(cns[cns > 1])))
  # Subsample high CNs to reduce computational overhead
  if(max_cns > 10){
    cns <- c(0:9, seq(from = 10, to = min(max_cns, 15), by = 5))
    if(max_cns > 49){
      cns <- c(cns, seq(from = 50, to = min(max_cns, 90), by = 10))
    }
    if(max_cns > 99){
      cns <- c(cns, seq(from = 100, to = max_cns, by = 100))
    }
  }else{cns <- 0:max_cns}
  positions_mafs <- depth(maxpeak = maxpeak, d = d_diff, P = test_ploidy, n = cns)
  depth_maf_posterior <- parametric_gmm_fit(data = depth_data, means = positions_mafs, variances = variances)
  depth_maf_responsibilities <- depth_maf_posterior/rowSums(depth_maf_posterior)
  depth_maf_responsibilities[is.na(rowSums(depth_maf_responsibilities)),] <- 0
  # Next we perform VAF fitting using the depth probabilities as a prior
  ## Create a list of vectors where each element represents a copy number
  ## and the elements are vectors of allelic states
  arr <- lapply(cns, testMAF, tp = test_tp)
  ## Fill vectors to same length
  max_len <- max(sapply(arr, length))
  arr <- lapply(arr, function(x){
    x <- c(x, rep(NA, times = max_len - length(x)))
  })
  ## Fit a mixture model to the VAF data for each assumption of copy number 
  ## where mixture components are allelic states
  arr <- lapply(arr, function(x){
    x <- parametric_gmm_fit(data = maf_df$values, means = x, variances = maf_variances)
    x <- data.table(x)
    # Take mean likelihoods for each bin
    x <- x[,lapply(.SD, mean), by = maf_df$ind]
    return(x)
  })
  
  # Next we use the depth likelihoods as a prior to scale the VAF likelihoods
  range_vals <- 1:nrow(depth_maf_responsibilities)
  for(i in 1:length(arr)){
    # Get fits for the ith copy number
    i_arr <- data.table::copy(arr[[i]])
    # Get range of depth values with and without an associated VAF value
    i_range <- range_vals[!range_vals %in% i_arr$maf_df]
    nomaf <- data.table("maf_df" = i_range)
    # Fill columns for VAF-less bins with NAs to match i_arr
    while(ncol(nomaf) < ncol(i_arr)){
      nomaf <- cbind.data.frame(nomaf, data.table(rep(NA, times = nrow(nomaf))))
      names(nomaf) <- names(i_arr)[1:length(names(nomaf))]
    }
    na_cols <- names(i_arr)[which(apply(arr[[i]], 2, function(x)all(is.na(x))))]
    i_arr <- rbind(i_arr, nomaf)
    i_arr <- data.frame(i_arr[order(maf_df)])
    if(length(na_cols)){
      i_arr[,-(which(names(i_arr) %in% na_cols))][i_range,-1] <- depth_maf_posterior[i_range,i]
    }else{
      i_arr[i_range,-1] <- depth_maf_posterior[i_range, i]
    }
    i_arr <- data.table(i_arr)
    arr[[i]] = (i_arr[,-1] * depth_maf_responsibilities[,i])
  }
  ## For each copy number assumption (element of arr), we return the 
  ## probability of the most likely allele state at that copy number for each
  ## observation, and then cbind it all together to get a matrix of likelihoods
  ## where each column is a copy number state
  joint_probs <- lapply(arr, function(x){apply(x, 1, max, na.rm = T)}) %>% do.call(cbind, .)
  
  joint_probs <- data.table(joint_probs)
  names(joint_probs) <- as.character(cns)
  maf_posteriors <- joint_probs %>% apply(1, max) %>% mean()
  maf_scores <- data.frame("maf" = mean(maf_posteriors), "tp" = test_tp, "ploidy" = test_ploidy, "lowest" = lowest)
  return(list("model" = maf_scores, "jp_tbl" = joint_probs))
}

maf_gmm_fit_subclonal <- function(depth_data, vaf_data, chr_vec, means, variances, maf_variances, maxpeak, ploidy = NA){
  ## Compute weights of initial fit
  proportions = compute_responsibilities(data = depth_data, means = means, variances = variances)
  proportions = colSums(proportions)/sum(colSums(proportions))
  
  ## Identify components with at least 1% of weight to use as boundaries for VAF fitting
  obs_cns <- as.numeric(names(proportions))
  max_cns <- ploidy + 3
  
  incl_pos <- means[obs_cns < max_cns]
  incl_weights <- proportions[obs_cns < max_cns]
  
  int_means <- means[obs_cns - round(obs_cns) == 0]
  
  
  #which.1pct <- proportions > 0.01
  #which.1pct <- names(proportions)[which.1pct]
  #which.1pct.min <- as.numeric(which.1pct) %>% min
  #which.1pct.max <- as.numeric(which.1pct) %>% max
  #onepct_pos <- means[names(means) %in% which.1pct.min:which.1pct.max]
  #onepct_weights <- proportions[names(proportions) %in% which.1pct.min:which.1pct.max]
  
  d_diff <- diff(int_means)[1]
  
  ## get peak fits
  peak_fits <- parametric_gmm_fit(data = depth_data, means = means, variances = variances) %>% apply(1, which.max)
  
  peak_fits <- names(means)[peak_fits]
  
  peak_in <- peak_fits %in% obs_cns
  
  ploidy_data <- data.frame("depth" = depth_data, "vaf" = vaf_data, CN = peak_fits, stringsAsFactors = F)
  
  lm_df <- data.frame("CN" = as.numeric(names(means)), "depth" = means)
  
  model = lm(CN ~ depth, lm_df)
  
  cov_char <- list(maxpeak = maxpeak, diff = d_diff)
  
  test_tp <- tp_diffploidy(cov_char, new_ploidy = ploidy)
  
  cov_char <- get_coverage_characteristics(test_tp, ploidy, maxpeak)
  
  maf_ind <- 1:nrow(ploidy_data)
  maf_ind <- maf_ind[!is.na(ploidy_data$vaf)]
  
  maf_list <- lapply(unmerge_mafs(as.character(ploidy_data$vaf[!is.na(ploidy_data$vaf)])), as.numeric)
  
  names(maf_list) <- maf_ind
  
  maf_df <- stack(maf_list)
  
  maf_df$ind <- as.numeric(levels(maf_df$ind)[maf_df$ind])
  
  cns <- obs_cns
  
  positions_mafs <- cns * d_diff + cov_char$homd
  
  depth_maf_posterior <- parametric_gmm_fit(data = depth_data, means = positions_mafs, variances = variances)
  
  depth_maf_responsibilities <- depth_maf_posterior/rowSums(depth_maf_posterior)
  
  depth_maf_responsibilities[is.na(rowSums(depth_maf_responsibilities)),] <- 0

  arr <- lapply(cns, testMAF_sc, tp = test_tp)
  
  print(arr)
  
  max_len <- max(sapply(arr, length))
  
  ## Fill vectors to same length
  arr <- lapply(arr, function(x){
    x <- c(x, rep(NA, times = max_len - length(x)))
    #names(x) <- 0:(max_len-1)
  })
  
  arr <- lapply(arr, function(x){
    #print(x)
    #x <- t(as.matrix(x))
    maf_variances <- match_kde_height(data = maf_df$values, means = x[!is.na(x)], sd = 0.03)
    x <- parametric_gmm_fit(data = maf_df$values, means = x, variances = maf_variances)
    #t <- mixpdf_function(means = arr[[4]], proportions = colSums(x/rowSums(x))/sum(colSums(x/rowSums(x))), sd = maf_variances)
    #data.frame(x = density(maf_df$values)$x, y = density(maf_df$values)$y, prob = t(density(maf_df$values)$x)$y) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_line(aes(x = x, y = prob, color = "GMM")) + geom_vline(xintercept = arr[[4]], linetype = 3) + xlab("VAF") + ylab("Density")
    #x = as.data.frame(x)
    #x$f <- maf_df$ind
    x <- data.table(x)
    x <- x[,lapply(.SD, mean), by = maf_df$ind]
    return(x)
  }
  )
  
  range_vals <- 1:nrow(depth_maf_responsibilities)
  for(i in 1:length(arr)){
    i_arr <- data.table::copy(arr[[i]])
    i_range <- range_vals[!range_vals %in% i_arr$maf_df]
    nomaf <- data.table("maf_df" = i_range)
    while(ncol(nomaf) < ncol(i_arr)){
      nomaf <- cbind.data.frame(nomaf, data.table(rep(NA, times = nrow(nomaf))))
      names(nomaf) <- names(i_arr)[1:length(names(nomaf))]
    }
    na_cols <- names(i_arr)[which(apply(arr[[i]], 2, function(x)all(is.na(x))))]
    i_arr <- rbind(i_arr, nomaf)
    i_arr <- data.frame(i_arr[order(maf_df)])
    if(length(na_cols)){
      i_arr[,-(which(names(i_arr) %in% na_cols))][i_range,-1] <- depth_maf_posterior[i_range,i]
    }else{
      i_arr[i_range,-1] <- depth_maf_posterior[i_range, i]
    }
    i_arr <- data.table(i_arr)
    arr[[i]] = (i_arr[,-1] * depth_maf_responsibilities[,i])
  }
  
  joint_probs <- lapply(arr, function(x){apply(x, 1, max, na.rm = T)}) %>% do.call(cbind, .)
  
  joint_probs <- data.table(joint_probs)
  names(joint_probs) <- as.character(cns)
  
  maf_posteriors <- joint_probs %>% apply(1, max) %>% mean()
  
  maf_scores <- data.frame("maf" = mean(maf_posteriors), "tp" = test_tp, "ploidy" = ploidy)
  
  return(list("model" = maf_scores, "jp_tbl" = joint_probs))
  
}

maf_gmm_fit_subclonal_prior <- function(depth_data, vaf_data, chr_vec, means, variances, maf_variances, maxpeak, ploidy = ploidy, tp = tp, cn_list){
  ## Compute weights of initial fit
  proportions = compute_responsibilities(data = depth_data, means = means, variances = variances)
  proportions = colSums(proportions)/sum(colSums(proportions))
  
  depth_vaf_df <- split(data.frame("d" = depth_data, "v" = vaf_data, "c" = chr_vec), f = chr_vec)
  
  curr_d_diff <- get_coverage_characteristics(tp, ploidy, maxpeak)$diff
  
  tot_cns <- as.numeric(names(table_vec(do.call(c, cn_list))))
  
  chr_joints <- list()
  for(chromosome in names(cn_list)){
    chr_data <- depth_vaf_df[[chromosome]]
    chr_means <- depth(maxpeak = maxpeak, d = curr_d_diff, P = ploidy, n = cn_list[[chromosome]])
    depth_maf_posterior <- parametric_gmm_fit(data = chr_data$d, means = chr_means, variances = variances)
    depth_maf_responsibilities <- depth_maf_posterior/rowSums(depth_maf_posterior)
    depth_maf_responsibilities[is.na(rowSums(depth_maf_responsibilities)),] <- 0
    cns <- as.numeric(names(chr_means))
    arr <- lapply(cns, testMAF_sc, tp = tp)
    maf_ind <- 1:nrow(chr_data)
    maf_ind <- maf_ind[!is.na(chr_data$v)]
    maf_list <- lapply(unmerge_mafs(as.character(chr_data$v[maf_ind])), as.numeric)
    names(maf_list) <- maf_ind
    maf_df <- stack(maf_list)
    maf_df$ind <- as.numeric(levels(maf_df$ind)[maf_df$ind])
    
    max_len <- max(sapply(arr, length))
    
    ## Fill vectors to same length
    arr <- lapply(arr, function(x){
      x <- c(x, rep(NA, times = max_len - length(x)))
      #names(x) <- 0:(max_len-1)
    })
    
    arr <- lapply(arr, function(x){
      #x <- t(as.matrix(x))
      #maf_variances <- match_kde_height(data = maf_df$values, means = x[!is.na(x)], sd = 0.03)
      
      x <- parametric_gmm_fit(data = maf_df$values, means = x, variances = maf_variances)
      #plot_density_gmm(maf_df$values, means = arr[[4]][1:3], sd = maf_variances, weights = colSums(x/rowSums(x, na.rm = T), na.rm = T)[1:3])
      
      #t <- mixpdf_function(means = arr[[4]], proportions = colSums(x/rowSums(x))/sum(colSums(x/rowSums(x))), sd = maf_variances)
      #data.frame(x = density(maf_df$values)$x, y = density(maf_df$values)$y, prob = t(density(maf_df$values)$x)$y) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_line(aes(x = x, y = prob, color = "GMM")) + geom_vline(xintercept = arr[[4]], linetype = 3) + xlab("VAF") + ylab("Density")
      #x = as.data.frame(x)
      #x$f <- maf_df$ind
      x <- data.table(x)
      x <- x[,lapply(.SD, mean), by = maf_df$ind]
      return(x)
    })
    
    range_vals <- 1:nrow(depth_maf_responsibilities)
    for(i in 1:length(arr)){
      i_arr <- data.table::copy(arr[[i]])
      i_range <- range_vals[!range_vals %in% i_arr$maf_df]
      nomaf <- data.table("maf_df" = i_range)
      while(ncol(nomaf) < ncol(i_arr)){
        nomaf <- cbind.data.frame(nomaf, data.table(rep(NA, times = nrow(nomaf))))
        names(nomaf) <- names(i_arr)[1:length(names(nomaf))]
      }
      na_cols <- names(i_arr)[which(apply(arr[[i]], 2, function(x)all(is.na(x))))]
      i_arr <- rbind(i_arr, nomaf)
      i_arr <- data.frame(i_arr[order(maf_df)])
      if(length(na_cols)){
        i_arr[,-(which(names(i_arr) %in% na_cols))][i_range,-1] <- depth_maf_posterior[i_range,i]
      }else{
        i_arr[i_range,-1] <- depth_maf_posterior[i_range, i]
      }
      i_arr <- data.table(i_arr)
      arr[[i]] = (i_arr[,-1] * depth_maf_responsibilities[,i])
    }
    
    joint_probs <- lapply(arr, function(x){apply(x, 1, max, na.rm = T)}) %>% do.call(cbind, .)
    
    colnames(joint_probs) <- cns
    
    to_add <- tot_cns[!as.character(tot_cns) %in% colnames(joint_probs)]
    
    to_add_dt <- data.table(matrix(0, nrow = nrow(joint_probs), ncol = length(to_add)))
    
    names(to_add_dt) <- as.character(to_add)
    
    joint_probs <- cbind(to_add_dt, joint_probs)
    
    setcolorder(joint_probs, as.character(tot_cns))
    
    joint_probs <- data.table(joint_probs)
    chr_joints <- c(chr_joints, list(joint_probs))
    names(chr_joints)[length(chr_joints)] <- chromosome
  }
  chr_joints <- rbindlist(chr_joints)
  
  maf_posteriors <- chr_joints %>% apply(1, max) %>% mean()
  
  maf_scores <- data.frame("maf" = mean(maf_posteriors), "tp" = tp, "ploidy" = ploidy)
  
  return(list("model" = maf_scores, "jp_tbl" = chr_joints))
  
}

na_or_true <- function(x){
  return((is.na(x)) | x)
}


plot_density_gmm <- function(data, means, weights, sd, ...){
  data <- data[data < max(means)]
  weights <- weights/sum(weights)
  den <- density(data, n = 2^16)
  pdf_fun <- mixpdf_function(means, weights, sd)
  den_pdf <- data.frame(x = den$x, y = den$y, prob = pdf_fun(den$x)$y)
  plot <- den_pdf %>% filter(x < max(means)) %>% ggplot(aes(x = x, y = y)) + geom_line() + geom_line(aes(x = x, y = prob, color = "Predicted")) + geom_vline(xintercept = means[means < max(data)], alpha = 0.1)
  return(plot)
  
}

lagged_df <- function(in_df){
  out <- rbindlist(list(in_df, in_df[1,]))
  out[nrow(out),] <- NA
  return(out[-1,])
}

staggered_seq <- function(to){
  out <- seq(from = 0, to = 10)
  if(to <= 10){
    return(out)
  }
  out <- c(out, seq(from = 15, to = min(max(15, ceiling(to/5)*5), 50), by = 5))
  if(to <= 50){
    return(out)
  }
  out <- c(out, seq(from = 60, to = min(max(60, ceiling(to/10)*10), 100), by = 10))
  if(to <= 100){
    return(out)
  }
  out <- c(out, seq(from = 150, to = min(max(150, ceiling(to/50)*50), 1000), by = 50))
  if(to <= 1000){
    return(out)
  }
  out <- c(out, seq(from = 1250, to = max(1250, ceiling(to/250)*250), by = 250))
  return(out)
}

n50_fun <- function(lengths){
  target <- sum(as.numeric(lengths))*0.5
  lengths <- sort(lengths)
  to_add <- lengths[1]
  s_sum = to_add
  while(s_sum < target){
    to_add <- lengths[1]
    lengths <- lengths[-1]
    s_sum <- s_sum + to_add
  }
  return(lengths[1])
}


len_fix <- function(x){
  if(length(x) == 0){
    NA
  }
  else{
    x
  }
}

weighted_median <- function(x, w){
  na_ind <- is.na(x)
  x <- x[!na_ind]
  w <- w[!na_ind]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  middling <- cumsum(w)/sum(w)
  which.mid <- max(which(middling <= 0.5))
  return(x[which.mid])
}
