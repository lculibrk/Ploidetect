# Used in ploidetect_cna_sc
get_good_var = function(mean1, mean2, base){
  sig = sqrt(((mean2 - mean1)^2)/(2 * log(base)))
  return(sig)
}

# Used in runiterativecompression_prob
compressdata_prob <- function(compress, criteria, dist_vec, subclones_discovered = F){
  ## Find length of input data
  reps <- lapply(compress, nrow)
  ## Get segment mappings
  ids = rep(x = 1:length(reps), times = unlist(reps))
  ## Record original raw distances
  dists_original <- dist_vec
  ## Get segmented distances
  dists <- aggregate(dist_vec, list(ids), FUN = mean)
  ## data.table of data to be segmented
  compress_original = rbindlist(compress)
  ## Average the data by segment
  compressed_compress <- compress_original[, lapply(.SD, mean), by = factor(ids)][,-1]
  ## Copy the averaged data so we can do stuff with it
  rel_liks <- data.table::copy(compressed_compress)
  ## Get maximum likelihood fit for each segment
  fits <- apply(rel_liks, 1, which.max)
  ## Look for regions which shift by at more than one component
  big_shift <- which(abs(diff(fits)) > 1)
  ## If we've already found subclones, compute transition likelihoods with adjacent segments included
  if(subclones_discovered){
    ## Get vector of states
    states <- apply(rel_liks, 1, which.max)
    ## Get maximum likelihood for each segment
    state_probs <- apply(rel_liks, 1, max)
    ## Get lagged data.table of data
    shifted <- lagged_df(rel_liks)
    ## Transform to data.frame
    t_liks <- as.data.frame(rel_liks)
    ## Vector of GMM fits
    state_vec <- 1:length(fits)
    ## Make vectors that return the current best & the following segment's best fit
    fit_val <- vapply(1:length(fits), function(x){t_liks[x,fits[x]]}, 0.01)
    fit_next <- vapply(1:length(fits), function(x){t_liks[x,c(fits, 1)[x+1]]}, 0.01)
    ## Bind them to make a lagged data.frame
    ## Does this state-by-state
    fit_df <- cbind(fit_val, fit_next)
    fit_df <- fit_df/rowSums(fit_df)
    fit_shifted <- vapply(1:length(states), function(x){t_liks[x + 1, states[x]]}, 0.01)
    fit_shifted_next <- vapply(1:length(states), function(x){t_liks[x+1,c(states, 1)[x+1]]}, 0.01)
    fit_shifted_df <- cbind(fit_shifted, fit_shifted_next)
    fit_shifted_df <- fit_shifted_df/rowSums(fit_shifted_df)
    ## Get differences in fits for transition "likelihoods"
    transition_liks <- rowSums(abs(fit_df - fit_shifted_df))
    transition_probs <- transition_liks[-length(transition_liks)]
    ## If NA, make them really big
    transition_probs[which(is.na(transition_probs))] <- dists$x[which(is.na(transition_probs))]
  }else{
    ## Vanilla transition likelihoods
    if(nrow(rel_liks) > 2){
      transition_probs <- rowSums(abs(apply(rel_liks, 2, diff)))
    }else if(nrow(rel_liks) == 2){
      ## This means there's only 2 segments here
      transition_probs <- sum(abs(apply(rel_liks, 2, diff)))
    }
  }
  ## TODO: Remove dependence on igraph, probably either trivial or a hair pulling exercise in self-hatred
  ## Initialize graph from data
  # Have to give edges a vector of numbers mapping edges to vertices, the code results in something like 
  # c(1, 2, 2, 3, 3, ... n-1, n-1, n) IIRC
  graph <- graph(edges = c(row.names(compressed_compress)[1],
                           rep(row.names(setDF(compressed_compress)[-c(1, nrow(compressed_compress)),]), each = 2),
                           row.names(setDF(compressed_compress))[nrow(setDF(compressed_compress))]), directed = F)
  ## Set edges to have the diffs attribute
  graph <- set_edge_attr(graph, name = "transition", value = transition_probs)
  ## Give vertices appropriate attributes
  graph <- set_vertex_attr(graph, name = "probs", value = compress)
  graph <- set_vertex_attr(graph, name = "npoints", value = unlist(reps))
  ## loop over all vertices
  sort_by_diff <- data.frame("vertex" = V(graph)$name)
  sort_by_diff$diff <- 0
  ## Apply transition likelihoods as edges
  edges <- edge_attr(graph, "transition", E(graph))
  ## We need to make a decision for which edge to delete for each vertex, since we always delete 
  ## at least one to ensure we only merge two adjacent bins, and not a long daisy-chain of them.
  e1 <- c(NA, edges)
  e2 <- c(edges, NA)
  ## Get magnitudes of transitions
  diff_vec <- abs(diff(dist_vec))
  ## For each point, we have the lagged and the led transitions
  d1 <- shift(dist_vec, type = c("lag"))
  d2 <- shift(dist_vec, type = c("lead"))
  ## Compare the two
  e_comp <- as.numeric(na_or_true(e1 > e2))
  ## Find ties, use the magnitude of the shifts to tiebreak
  e_eq <- which(e1 == e2)
  d_comp <- as.numeric(na_or_true(d1 > d2))
  e_comp[e_eq] <- d_comp[e_eq]
  ## This checks if we're on iteration 1 of the recursive segmentation
  if(all(unlist(reps) == 1)){
    ## I'll have to run this line-by-line to remember exactly what this does
    ## But I think this checks for remaining sequentially linked vertices that weren't found
    ## in the previous step
    preserve <- unique(1:length(V(graph)) - (1-e_comp))
    sequentials <- preserve[which(na_or_true((preserve - shift(preserve) == 1)))]
    sequentials <- sequentials[!sequentials %in% c(1, length(edges) + 1)]
    preserve <- preserve[!preserve %in% (sequentials - e_comp[sequentials])]
    to_del <- V(graph)[!V(graph) %in% preserve]
    graph <- delete_edges(graph, to_del)
  }else{
    # Get incident vertices of cases where both edges are beyond the threshold
    both_ind <- na_or_true(e1 > criteria) & na_or_true(e2 > criteria)
    both_ind <- which(both_ind & (unlist(reps) > 1))
    both <- unique(sort(pmax(1, c(both_ind, both_ind - 1))))
    both <- both[both < length(e1)]
    ## Check for vertices where there is only one bin in the segment
    lenone <- which(unlist(reps) == 1)
    big_shift <- big_shift[!big_shift %in% lenone]
    ## Preserve the edge that has the lower differential probability
    preserve <- unique(lenone - (1 - e_comp[lenone]))
    ## As above, i think this checks for sequentially linked vertices
    sequentials <- preserve[which(na_or_true((preserve - shift(preserve) == 1) & (-(preserve - shift(preserve, type = "lead")) == 1)))]
    sequentials <- sequentials[!sequentials %in% c(1, length(edges) + 1)]
    sequentials <- sequentials[sequentials > 0]
    preserve <- preserve[!preserve %in% (sequentials - na_or_true(edges[sequentials - 1] > edges[sequentials + 1]))]
    
    # e1 is left edge, e2 is right edge for each vertex
    # Here we subtract the vertex IDs by an integer-boolean of whether e1 > e2
    # If e1 is the larger edge, then this evaluates to TRUE, so we get the value vetex_id - 1,
    # which when translated to edge IDs, is the left edge, resulting in the left edge being deleted.
    single <- as.numeric(na.omit(unique((1:length(e1)) - e_comp)))
    single <- single[single > 0]
    single <- unique(c(single, which(edges > criteria)))
    to_delete <- unique(c(both, single, big_shift))
    to_delete <- to_delete[!to_delete %in% preserve]
    ## Delete edges we don't want contracted
    graph <- delete_edges(graph, to_delete)
  }
  
  # Get list of all vertex pairs to merge
  tomerge <- ends(graph, E(graph))
  merging <- as.numeric(c(t(tomerge)))
  vertices <- components(graph)$membership
  # Merge connected vertices
  graph <- contract.vertices(graph, mapping = vertices, vertex.attr.comb = list("probs" = function(...){
    elements = c(...)
    return(rbindlist(elements))
  }, "npoints" = "sum", "name" = "first"))
  
  ## Reconstruct the data.frame we began with
  out_probs = get.vertex.attribute(graph, "probs")
  n = get.vertex.attribute(graph, "npoints")
  ## Get compressed probs
  compress <- out_probs
  i_segs <- 1:length(compress)
  i_segs <- rep(i_segs, times = unlist(lapply(compress, nrow)))
  ## Verbosity message
  if(verbose){
    message(paste0("Iteration complete, segment count = ", length(compress)))
  }
  return(compress)
}

# Used in ploidetect_prob_segmentator
runiterativecompression_prob <- function(data, segmentation_threshold = segmentation_threshold, verbose = F, subclones_discovered = F){
  ## Verbosity statement
  if(verbose){
    message("Running iterative compression to segment read depth data")
  }
  ## Define segments and the raw distance
  segments <- data$seg
  dist = data$dist
  
  ## Separate data from metadata
  data <- data[,-which(names(data) %in% c("chr", "seg", "dist"))]
  ## Verbosity statement
  if(verbose){
    message(paste0("Initial segment count: ", length(unique(segments))))
  }
  ## Set initial conditions - i.e. not converged, generate initial data.table, and split by segment.
  converged <- F
  compress <- data.table(data)
  compress <- split(compress, f = segments)
  ## If a single segment was input
  if(length(compress) == 1){
    converged = T
  }
  
  while(!converged){
    ## Record number of segments pre-segmentation
    windows <- length(compress)
    ## Run an iteration of compression segmentation using the input threshold
    compress <- compressdata_prob(compress, criteria = segmentation_threshold, dist_vec = dist, subclones_discovered = subclones_discovered)
    ## Check if we've converged or if there's only one segment left
    if(length(compress) == windows | length(compress) == 1){
      converged <- T
    }
  }
  ## Output the segment mappings
  segs <- rep(1:length(compress), times = lapply(compress, nrow))
  return(segs)
}

# Used in segment_subclones, ploidetect_cna_sc
ploidetect_prob_segmentator <- function(prob_mat, ploidy, chr_vec, seg_vec, dist_vec, verbose = T, lik_shift = 1, subclones_discovered = F){
  ## Verbosity statement
  if(verbose){
    message("Performing segmentation and copy number variant calling")
  }
  ## Standardizing the input data
  if(!all(c("chr", "seg") %in% names(prob_mat))){
    prob_mat <- cbind(prob_mat, "chr" = chr_vec, "seg" = seg_vec)
  }
  
  ## Keep a copy of the initial matrix
  orig_mat <- data.table::copy(prob_mat)
  
  ## Make a copy that will be converted to posteriors
  resp_mat <- data.table::copy(prob_mat)
  
  ## remove the chr, seg columns from resp_mat to calc the posteriors
  resp_mat <- resp_mat[,-((ncol(resp_mat)-1):ncol(resp_mat))]
  
  ## Calculate posteriors
  resp_mat <- resp_mat/rowSums(resp_mat)
  
  ## Correct divisions by zero
  resp_mat[which(is.na(resp_mat[,1])),] <- 0
  
  ## Take the mean over the input segments
  compressed_mat <- resp_mat[, lapply(.SD, mean), by = list(orig_mat$chr, orig_mat$seg)]
  
  ## Remove chr, seg columns that were introduced in the previous line
  compressed_mat <- compressed_mat[,-c(1:2)]
  
  ## Checks if the matrix has shrunk at all i.e. if the input was previously segmented
  if(nrow(compressed_mat) < nrow(orig_mat)){
    ## get the MLE GMM component for each segment
    states <- apply(compressed_mat, 1, which.max)
    ## Create a leading/lagging data.table to allow for checking for transitions between components
    state_transitions <- data.frame("leading" = states[-1], "lagging" = states[1:(length(states)-1)])
    ## Split the data.table at each row
    state_transitions <- split(state_transitions, f = 1:nrow(state_transitions))
    ## Convert each resultant list element to a vector
    state_transitions <- lapply(state_transitions, unlist)
    ## Get diff of each column in each row
    all_transitions <- apply(abs(apply(compressed_mat, 2, diff)), 1, sum)
    ## Get transition threshold - either 50th percentile of possible transitions, or 1.5, whichever is lower.
    transition_threshold <- min(1.5, quantile(all_transitions[all_transitions > 0.5], prob = 0.5))
    ## If the data is uncompressed, use the shift specified in the function call
  }else{transition_threshold = lik_shift}
  ## Build data.table for input to the segmentation function
  prob_mat$chr <- chr_vec
  prob_mat$seg <- seg_vec
  prob_mat$dist <- dist_vec
  ## Set NAs to 0
  prob_mat[which(is.na(prob_mat[,1])),] <- 0
  ## Split by chromosome
  datsplit <- split(as.data.frame(prob_mat), prob_mat$chr)
  ## Verbosity statement
  if(verbose){
    message("Performing segmentation of copy number data")
  }
  ## Run runiterativecompression_prob on data
  compressedalldat <- unlist(lapply(datsplit, runiterativecompression_prob, segmentation_threshold = transition_threshold, verbose = verbose, subclones_discovered = subclones_discovered))
  return(compressedalldat)
}
# Used in segment_subclones
# Misnamed - used to use BAF & depth (hence joint), but turns out depth only is more accurate. 
subclonal_joint_probs <- function(new_seg_data, tp, ploidy, maxpeak, clonal_positions, vaf_variance, fracs = c(0.25, 0.5), in_sd = NA){
  ## Unpack mapping of CN values to read depths
  predictedpositions <- get_coverage_characteristics(tp, ploidy, maxpeak)$cn_by_depth
  ## Generate an LM for predicting subclonal CNVs from depth
  train_df <- data.frame("segment_depth" = clonal_positions, "CN" = as.numeric(names(predictedpositions)))
  model <- lm(CN ~ segment_depth, train_df)
  ## Get depth differential from TP/ploidy estimate
  d_diff <- get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)$diff
  ## Pick which fraction we will use for subclonal CNA calling
  ## Could be faster by not concatenating to vectors/lists, but this doesn't loop many times
  int_liks <- c()
  variances <- c()
  positions_vectors <- list()
  for(subcl_frac in fracs){
    ## Call subclonal CNVs by position as an estimate
    new_seg_data$CN <- round(predict(model, new_seg_data)/(10*subcl_frac), digits = 1)*10*subcl_frac
    ## Limit subclonal CNV calling to ploidy+3 maximum, not really meaningful past that point.
    subcl_pos <- seq(from = min(round(new_seg_data$CN)), to = max(round(new_seg_data$CN)), by = subcl_frac)
    ## Generate a vector of depths
    subcl_pos <- structure(depth(maxpeak, d = d_diff, P = ploidy, n = seq(from = min(round(new_seg_data$CN)), to = max(round(new_seg_data$CN)), by = subcl_frac)), names = subcl_pos)
    ## Make sure the positions include CN from 0 to 10
    to_add <- (0:10)[which(!0:10 %in% names(subcl_pos))]
    subcl_pos <- sort(c(subcl_pos, depth(maxpeak, d = d_diff, P = ploidy, n = to_add)))
    ## Filter for 0 to ploidy + 3
    subcl_pos <- subcl_pos[!(as.numeric(names(subcl_pos)) > ploidy + 3 & as.numeric(names(subcl_pos)) - floor(as.numeric(names(subcl_pos))) != 0)]
    ## Store positions
    positions_vectors <- c(positions_vectors, list(subcl_pos))
    ## Variance to be used can be given (I don't think it ever is, but I won't test removing this block yet)
    if(is.na(in_sd)){
      ## Variance estimation
      ## Compute variance by segment, grab only largest 50% of segments or top ten, whichever is smaller
      segged_sd = new_seg_data[CN - round(CN) == 0]
      segged_sd[,n:=.N, by = list(chr, segment)]
      segged_sd = mean(segged_sd[,.(var = sd(corrected_depth), n = first(n)), by = list(chr, segment)][order(n, decreasing = T)][1:min(10, .N/2)]$var)
      ## Use the previous estimate to get variance by comparison to a KDE
      segged_sd <- match_kde_height(new_seg_data$segment_depth, means = subcl_pos, sd = segged_sd)
    }else{segged_sd <- in_sd}
    ## Store sd
    variances <- c(variances, segged_sd)
    subcl_positions <- as.numeric(names(subcl_pos)) - round(as.numeric(names(subcl_pos))) != 0
    segged_sd <- rep(segged_sd, length.out = length(subcl_pos))
    segged_sd[subcl_positions] <- segged_sd[subcl_positions]/2
    
    ## Fit segmented depth to new subclonal-aware GMM
    subcl_fit <- parametric_gmm_fit(new_seg_data$segment_depth, means = subcl_pos, variances = segged_sd)
    
    ## Get responsibilities
    subcl_resps <- subcl_fit/rowSums(subcl_fit)
    subcl_resps[apply(subcl_resps, 1, function(x)all(is.na(x))),] <- 0
    
    ## Compute likelihood values
    new_max_resps <- apply(subcl_resps * subcl_fit, 1, sum)
    
    ## Get mean likelihood
    int_lik <- mean(new_max_resps[new_seg_data$CN - round(new_seg_data$CN) != 0], na.rm = T)
    
    ## Add to list
    int_liks <- c(int_liks, int_lik)
  }
  ## Pick positions and SD based on which model gave the better likelihoods
  subcl_pos <- positions_vectors[[which.max(int_liks)]]
  segged_sd <- variances[[which.max(int_liks)]]
  ## Fit a new gmm using the positions
  subclonal_probs <- data.table(parametric_gmm_fit(data = new_seg_data$segment_depth, means = subcl_pos, variances = segged_sd))
  colnames(subclonal_probs) <- names(subcl_pos)
  ## Bit of holdover code to put it into the list
  subclonal_probs <- list("jp_tbl" = subclonal_probs)
  return(list("fraction" = fracs[which.max(int_liks)], "probs" = subclonal_probs, "segged_sd" = segged_sd))
}

# Used in ploidetect_cna_sc
segment_subclones <- function(new_seg_data, predictedpositions, depth_variance, vaf_variance, maxpeak, tp, ploidy){
  ## Generate an LM for predicting subclonal CNVs from depth
  train_df <- data.frame("segment_depth" = predictedpositions, "CN" = as.numeric(names(predictedpositions)))
  model <- lm(CN ~ segment_depth, train_df)
  ## Compute joint probabilities
  subclonal_probs <- subclonal_joint_probs(new_seg_data = new_seg_data, tp = tp, ploidy = ploidy, maxpeak = maxpeak, clonal_positions = predictedpositions, vaf_variance = 0.03)
  ##Store fraction
  subclonal_fraction <- subclonal_probs$fraction
  subclonal_variance <- subclonal_probs$segged_sd
  subclonal_probs <- subclonal_probs$probs
  ## Generate segments from subclonal probability matrix using segmentation by compression
  newvelle_segs <- ploidetect_prob_segmentator(prob_mat = subclonal_probs$jp_tbl, ploidy = ploidy, chr_vec = new_seg_data$chr, seg_vec = unlist(lapply(split(1:nrow(new_seg_data), new_seg_data$chr), function(x)1:length(x))), dist_vec = new_seg_data$corrected_depth, subclones_discovered = T, lik_shift = 1.5)
  new_seg_data$segment <- newvelle_segs
  new_seg_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  ## Map copy number states to segments from GMM
  new_seg_data$CN <- as.numeric(names(subclonal_probs$jp_tbl)[apply(subclonal_probs$jp_tbl, 1, which.max)])
  ## Get CN states not found in the GMM & round to the decided-on subclonal fraction used
  new_seg_data[is.na(CN)]$CN <- round(predict(model, data.frame("segment_depth" = new_seg_data[is.na(CN)]$segment_depth))/(subclonal_fraction*10), 1)*(subclonal_fraction*10)
  ## Get catalogue of all CN states in each chromosome
  CNA_list <- new_seg_data[,c("chr", "CN")]
  CNA_list <- unique(CNA_list)
  ## Remove subclonal CNVs above ploidy+3
  CNA_list$CN[CNA_list$CN > ploidy + 3] <- round(CNA_list$CN[CNA_list$CN > ploidy + 3])
  ## Get filtered predicted positions for each CN
  comp_pos <- depth(maxpeak = maxpeak, d = get_coverage_characteristics(tp, ploidy, maxpeak)$diff, P = ploidy, sort(unique(CNA_list$CN)))
  ## Remove negative CNs that might be included due to noise
  CNA_list$CN <- pmax(0, CNA_list$CN)
  ## Remove duplicate copy number states from each chromosome
  CNA_list <- split(CNA_list$CN, f = CNA_list$chr)
  CNA_list <- lapply(CNA_list, function(x)unique(sort(x)))
  ## Re-calculate likelihoods based on filtered positions
  new_jp_tbl <- list("jp_tbl" = data.table(parametric_gmm_fit(data = new_seg_data$segment_depth, means = comp_pos, variances = depth_variance)))
  #### Generate GMM fits for subclonal copy number calls from joint probability matrix
  segmented_subclonal_probs <- new_jp_tbl$jp_tbl
  ## Add grouping variables to joint probabilities
  segmented_subclonal_probs$chr = new_seg_data$chr
  segmented_subclonal_probs$segment = newvelle_segs
  ## Summarise probabilities
  segmented_subclonal_probs <- segmented_subclonal_probs[, lapply(.SD, mean), by = list(chr, segment)]
  ## Pull out state vector
  states <- names(segmented_subclonal_probs)[-(1:2)]
  ## Pull out grouping variables
  seg_info <- segmented_subclonal_probs[,1:2]
  segmented_subclonal_probs <- segmented_subclonal_probs[,-(1:2)]
  ## Get maximum likelihood state for each segment
  segmented_subclonal_probs[, call := states[which.max(.SD)], by = 1:nrow(segmented_subclonal_probs)]
  ## Add back segment variables and select only those + calls
  segmented_subclonal_probs <- cbind(seg_info, segmented_subclonal_probs)
  segmented_subclonal_probs <- segmented_subclonal_probs[,c("chr", "segment", "call")]
  ## Left join state calls with data
  new_seg_data <- data.table(new_seg_data)
  new_seg_data[segmented_subclonal_probs, on = list(chr, segment), call := as.numeric(call)]
  new_seg_data$CN <- new_seg_data$call
  ## Get "parent" CNAs
  new_seg_data$parent_cns <- round(new_seg_data$CN/5, digits = 1)*5
  ## Compute size of subclonal events called
  putative_subclones <- new_seg_data[,.(size = last(end) - first(pos), CN = first(CN), first.n = first(.I), last.n = last(.I)), by = list(chr, segment)]
  ## Filter for subclonal events
  subcl_vec <- which(putative_subclones$CN != round(putative_subclones$CN))
  ## Filter for "good" subclones - Heuristic
  ## Basically make sure that we don't have runs of only subclonal CNVs
  blacklist_subcl <- putative_subclones[subcl_vec[na_or_true(shift(putative_subclones$CN, type = "lag")[subcl_vec] == shift(putative_subclones$CN, type = "lead")[subcl_vec] & shift(putative_subclones$chr, type = "lag")[subcl_vec] == shift(putative_subclones$chr, type = "lead")[subcl_vec])]]
  same_cn_preceding <- putative_subclones[abs(CN - shift(CN, type = "lead")) < 1 & na_or_true(chr != shift(chr, type = "lag"))][CN != round(CN)]
  same_cn_following <- putative_subclones[abs(CN - shift(CN, type = "lag")) < 1 & na_or_true(chr != shift(chr, type = "lead"))][CN != round(CN)]
  ## Remove subclones which comprise fewer than 5Mb of sequence - can probably pump this up a bit
  blacklist_subcl <- unique(rbind(blacklist_subcl, same_cn_following, same_cn_preceding)[size < 5e+06])
  blacklist_verts <- sort(c(blacklist_subcl$first.n - 1, blacklist_subcl$last.n))
  ## Fix subclonal segments that might be noisy using igraph
  edge_vec <- c(1, rep(2:(nrow(new_seg_data)-1), each = 2), nrow(new_seg_data))
  g <- graph(edges = edge_vec, directed = F)
  g <- set_vertex_attr(g, name = "chr", value = new_seg_data$chr)
  todel <- c(which(abs(diff(new_seg_data$CN)) > 0), 
             which(!na_or_true(shift(new_seg_data$chr, type = "lead") == new_seg_data$chr))
  )
  todel <- todel[!todel %in% blacklist_verts]
  g <- delete_edges(g, todel)
  
  ## Apply fixed segments to data
  new_seg_data$segment <- components(g)$membership
  new_seg_data[,segment_depth:=median(corrected_depth), by = list(chr, segment)]
  new_seg_data[,CN:=as.numeric(names(comp_pos)[apply(parametric_gmm_fit(data = segment_depth, means = comp_pos, variances = depth_variance), 1, which.max)])]
  return(list("data" = list(new_seg_data), "fraction" = subclonal_fraction, "subclonal_variance" = subclonal_variance))
}

## Main segmentation function
#' @export
ploidetect_cna_sc <- function(all_data, segmented_data, tp, ploidy, maxpeak, verbose = T, min_size = 1, simp_size = 100000, max_iters = Inf){
  ## Get estimated differential depth
  d_diff <- get_coverage_characteristics(tp, ploidy, maxpeak)$diff
  ## Get estimated positions for up to 1000-fold amplification
  predictedpositions <- depth(maxpeak = maxpeak, d = d_diff, P = ploidy, n = 0:1000)
  ## Filter positions for ones which actually exist
  predictedpositions <- predictedpositions[predictedpositions < max(segmented_data$corrected_depth)]
  ## Ensure predictedpositions has at least CNs 0-10
  if(any(!as.character(0:10) %in% names(predictedpositions))){
    missing = (0:10)[!as.character(0:10) %in% names(predictedpositions)]
    predictedpositions = sort(c(predictedpositions, depth(maxpeak, d_diff, ploidy, n = missing)))
  }
  ## Correct for X chromosome being single-copy in males that was otherwise normalized out during preprocessing
  ## Check if the case is a male
  sizes = all_data$end - all_data$pos
  autosomes = sizes[!all_data$chr %in% c("X", "Y")]
  sex_chrs = sizes[all_data$chr %in% c("X", "Y")]
  expected_x = median(autosomes)/2
  sex_fit = parametric_gmm_fit(sex_chrs, c(expected_x, median(autosomes)), get_good_var(expected_x, median(autosomes), base = 5))
  if(which.max(colSums(sex_fit)) == 1){
    ## Male
    all_data$tumour[all_data$chr == "X"] = all_data$tumour[all_data$chr == "X"]/2
    segmented_data[chr == "X"]$corrected_depth = segmented_data[chr == "X"]$corrected_depth/2
  }
  
  ## Estimate variance based on KDE matching
  variance <- density(segmented_data$corrected_depth)$bw
  variance <- match_kde_height(segmented_data$corrected_depth, means = predictedpositions, sd = variance)
  
  ##Compute KDE weights
  proportions <- compute_responsibilities(segmented_data$corrected_depth, means = predictedpositions, variances = variance)
  proportions <- colSums(proportions)/sum(colSums(proportions))
  
  ## Recompute ploidy & point of highest density
  ploidy <- as.numeric(names(proportions)[which.max(proportions)])
  maxpeak <- predictedpositions[which.max(proportions)]
  
  ## Recompute var based on parameters
  variance <- match_kde_height(segmented_data$corrected_depth, means = predictedpositions, sd = variance, comparison_point = maxpeak)
  
  ## Compute likelihoods based on GMM fit
  joint_probs <- list("jp_tbl" = data.table(parametric_gmm_fit(segmented_data$corrected_depth, means = predictedpositions, variances = variance)))
  
  ## Get the 80th percentile of likelihood shifts as a very low bar for similarity
  lik_shift <- quantile(rowSums(abs(apply(joint_probs$jp_tbl, 2, diff))), prob = 0.8)
  
  ## Segment genome based on 80th percentile of likelihoods
  clonal_cnas <- ploidetect_prob_segmentator(prob_mat = joint_probs$jp_tbl, ploidy = ploidy, chr_vec = segmented_data$chr, seg_vec = unlist(lapply(split(1:nrow(segmented_data), segmented_data$chr), function(x)1:length(x))), dist_vec = segmented_data$corrected_depth, lik_shift = lik_shift)
  clonal_cna_data <- segmented_data
  clonal_cna_data$segment <- clonal_cnas
  
  
  ## Convert to data.table for efficiency
  clonal_cna_data <- data.table(clonal_cna_data)
  
  ## Compute segmented depth
  clonal_cna_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  
  ## Call integer copy numbers
  clonal_cna_data[,CN:=round(cn_from_dp(segment_depth, maxpeak, tp, ploidy))]
  
  ## Subclonal-aware segmentation
  subcl_seg <- segment_subclones(new_seg_data = clonal_cna_data, predictedpositions = predictedpositions[1:11], depth_variance = variance, vaf_variance = 0.06, maxpeak = maxpeak, tp = tp, ploidy = ploidy)
  ## Extract estimated fraction for subclones; 0.25 or 0.5
  subclonal_fraction <- subcl_seg$fraction
  ## Extract SD used for previous segmentation step
  subclonal_variance <- subcl_seg$subclonal_variance
  ## Extract segment data.table
  subcl_seg <- subcl_seg$data[[1]]
  ## Get which segments are predicted to be subclonal
  subcl_seg$subclonal <- (as.numeric(subcl_seg$CN) - round(subcl_seg$CN)) != 0
  ## Get segment mappings
  seg_mappings = subcl_seg[,.(pos = first(pos), end = last(end), segment_depth = first(segment_depth), call = first(call), n = .N), by = list(chr, segment)]
  ## Get clonal positions & use regression to get fractional positions
  predictedpositions = depth(maxpeak = maxpeak, d = d_diff, P = ploidy, n = 0:10)
  cn_df = data.frame(cn = 0:10, segment_depth = predictedpositions)
  cn_fit = lm(cn ~ segment_depth, data = cn_df)
  ## Get fractional CN calls
  seg_mappings$fine_call = round(predict(cn_fit, seg_mappings), 2)
  ## Segment mappings by chromosome
  ## occasionally local variation or extremely low abundance subclonal CNV may result in segments
  ## with an apparent fractional copy number slightly above or below their integer value. Here we 
  ## Record the fractional CNs for each segment to "nudge" the means of GMM fits on a per-segment
  ## basis later on, so that these slight shifts don't cause inflated transitions
  chr_mappings <- split(seg_mappings, f = seg_mappings$chr)
  individual_pos <- lapply(chr_mappings, function(x){
    sort(na.omit(unique(c(names(predictedpositions), unique(x$call)))))
  })
  pos_list <- lapply(1:nrow(seg_mappings), function(x){
    wp <- individual_pos[[seg_mappings$chr[x]]]
    wp[which(wp == seg_mappings$call[x])] <- seg_mappings$fine_call[x]
    as.numeric(wp)
  })
  
  ## Segment according to the position list obtained directly above
  refined_liks <- maf_gmm_fit_subclonal_prior_segments(depth_data = subcl_seg$segment_depth, vaf_data = subcl_seg$maf, chr_vec = subcl_seg$chr, means = predictedpositions, variances = variance, maf_variances = 0.06, maxpeak = maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = seg_mappings)$jp_tbl
  subcl_seg$segment <- ploidetect_prob_segmentator(prob_mat = refined_liks, ploidy = ploidy, chr_vec = subcl_seg$chr, seg_vec = subcl_seg$segment, dist_vec = subcl_seg$segment_depth, lik_shift = 1.5, subclones_discovered = T)
  subcl_seg[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  seg_mappings = subcl_seg[,.(pos = first(pos), end = last(end), segment_depth = first(segment_depth), call = first(call), n = .N), by = list(chr, segment)]
  
  ## Re-compute the positions to get a more accurate view of the subclones present and re-fit the nudged GMM
  cn_df = data.frame(cn = 0:10, segment_depth = predictedpositions)
  cn_fit = lm(cn ~ segment_depth, data = cn_df)
  seg_mappings$fine_call = round(predict(cn_fit, seg_mappings), 2)
  chr_mappings <- split(seg_mappings, f = seg_mappings$chr)
  individual_pos <- lapply(chr_mappings, function(x){
    sort(na.omit(unique(c(names(predictedpositions), unique(x$call)))))
  })
  pos_list <- lapply(1:nrow(seg_mappings), function(x){
    wp <- individual_pos[[seg_mappings$chr[x]]]
    wp[which(wp == seg_mappings$call[x])] <- seg_mappings$fine_call[x]
    as.numeric(wp)
  })
  refined_liks <- maf_gmm_fit_subclonal_prior_segments(depth_data = subcl_seg$segment_depth, vaf_data = subcl_seg$maf, chr_vec = subcl_seg$chr, means = predictedpositions, variances = variance, maf_variances = 0.06, maxpeak = maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = seg_mappings)$jp_tbl
  ## Output calls
  subcl_seg$call = apply(refined_liks, 1, function(x)names(refined_liks)[which.max(x)])
  seg_mappings <- split(seg_mappings, seg_mappings$chr)
  
  ## Get most common CN to estimate BAF variance
  ## TODO: Check the maf_gmm_fit_subclonal_prior_segments function to see if this is even needed
  common_call <- names(which.max(table_vec(subcl_seg$CN)))
  common_maf_means <- testMAF(as.numeric(common_call), tp)
  maf_variance <- match_kde_height(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, sd = 0.03)
  
  ## Subclonal positions
  obs_pos <- as.numeric(names(table_vec(subcl_seg$CN)))
  obs_pos <- sort(c(obs_pos, (0:10)[!0:10 %in% obs_pos]))
  subcl_seg$CN[subcl_seg$CN < 0] <- 0
  subcl_cn <- as.numeric(names(table_vec(subcl_seg$CN)))
  subcl_pos <- depth(maxpeak, d_diff, ploidy, subcl_seg$CN)
  
  ## Flag subclones
  subcl_seg$subclonal <- (as.numeric(subcl_seg$CN) - round(subcl_seg$CN)) != 0
  
  ## Re-compute BAF variance
  ## TODO: test if this actually does anything on top of doing this a few lines up
  common_call <- names(which.max(table_vec(subcl_seg$CN)))
  common_maf_means <- testMAF(as.numeric(common_call), tp)
  maf_variance <- match_kde_height(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, sd = 0.03)
  
  ## Get subclone calling LM
  ## TODO: test if this actually does anything different compared to earlier model
  split_segs <- split(subcl_seg, f = subcl_seg$chr)
  model <- lm(CN ~ corrected_depth, data = data.frame(CN = as.numeric(names(predictedpositions)), corrected_depth = predictedpositions))
  
  ## Get chromosome-wise seen positions
  individual_pos <- lapply(split_segs, function(x){
    vec <- as.numeric(names(table_vec(x$CN)))
    max_val <- max(ceiling(vec))
    ind_pos <- sort(c(vec, (0:max_val)[!0:max_val %in% vec]))
    cns <- names(table_vec(round(predict(model, x))))
    ind_pos <- sort(as.numeric(unique(c(ind_pos, cns))))
    ind_pos[ind_pos < 0] <- 0
    return(unique(ind_pos))
  })
  
  ## Get the number of bins merged into each meta-bin used in TC/Ploidy calling
  unaltered <- ploidetect_preprocess(all_data, simplify = T, simplify_size = simp_size, verbose = F)$merged
  ## Get un-merged bins with mappability/gc-bias normalization applied
  all_data_preprocessed <- ploidetect_preprocess(all_data, simplify = F, simplify_size = NA, verbose = T)
  all_data_preprocessed <- all_data_preprocessed$x
  ## Get segment mappings 
  initial_segment_mappings <- subcl_seg[,.(pos = first(pos), end = last(end), CN = first(CN)), by = list(chr, segment)]
  ## Map segments to initial data
  reduced_mappings <- initial_segment_mappings[,-"end"][all_data_preprocessed, on = c("chr", "pos"), roll = Inf]
  ## Get segment depth
  reduced_mappings[, segment_depth := median(corrected_depth), by = list(chr, segment)]
  ## Get the bins-to-metabins merging numbers for the coarse-to-fine part of segmentation
  iterations <- round(unaltered/2^(1:ceiling(log2(unaltered))))
  ## Ensure that the final set of coarse-to-fine iterations includes one iteration of 2 and one of 1 (raw)
  if(!all(c(1, 2) %in% iterations)){
    iterations <- c(iterations[1:min(which(iterations == 1) - 1)], 2, 1)
  }
  ## Re-estimate maxpeak & ploidy after the shifting around of the data that's been done
  ## Ploidy can sometimes change from the initial estimate, which is why we do this
  ploidy <- as.numeric(names(subcl_pos)[which.min(abs(density(subcl_seg$segment_depth, n = 2^16)$x[which.max(density(subcl_seg$segment_depth, n = 2^16)$y)] - subcl_pos))])
  maxpeak <- density(subcl_seg$segment_depth, n = 2^16)$x[which.max(density(subcl_seg$segment_depth, n = 2^16)$y)]
  cn_positions = get_coverage_characteristics(tp, ploidy, maxpeak)$cn_by_depth
  
  ## To estimate the true "peak" of read depth we filter down for regions that are half a predicted copy number
  ## away from the initially estimated maxpeak value
  closeness <- abs(subcl_seg$segment_depth - maxpeak)
  maxpeak_segments <- unique(subcl_seg[which(closeness < diff(predictedpositions)[1]/2), c("chr", "segment")])
  maxpeak_segments$mp <- T
  
  ## Re-compute subclonal positions
  subcl_pos <- depth(maxpeak = maxpeak, d = get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)$diff, P = ploidy, n = obs_pos)
  ## Get estimate of maxpeak for the non-merged data
  maxpeak_base <- maxpeak/max(1, (unaltered - 1))
  ## Get expected coverage characteristics for non-merged data
  base_characteristics <- get_coverage_characteristics(tp, ploidy, maxpeak_base)
  ## Initialize data.tables for the coarse-to-fine, since sometimes it aborts early
  previous_segment_mappings <- data.table::copy(initial_segment_mappings)
  overseg_mappings = data.table::copy(subcl_seg)
  overseg_mappings = overseg_mappings[,.(chr = chr, segment = segment, pos = pos, CN = CN, merge_vec = 1:.N, end = end, corrected_depth = corrected_depth, maf = maf, gc = 0.5, segment_depth = segment_depth, fine_call = CN, flagged = F, match = F, call = CN)]
  seg_lens <- previous_segment_mappings[,.("pos" = first(pos), "end" = last(end)), by = list(chr, segment)][,.(diff=end-pos)]$diff
  seg_lens <- seg_lens[seg_lens > 0]
  current_n50 <- n50_fun(seg_lens)
  current_median_length <- median(seg_lens)
  subclonal_seg_mappings <- setnames(rbindlist(seg_mappings), old = "call", new = "CN")
  ## Begin coarse-to-fine segmentation

  ## Set first iteration
  i = 1

  ## Exit condition
  condition = T
  while(condition){
    ## Record the "mergings" done to the input data in this iteration
    val <- iterations[i]
    ## If we've exceeded the maximum specified iterations, exit and return segments
    if(i > max_iters){
      condition = F
      break
    }
    ## If we've exceeded the maximum resolution, exit and return segments
    if(val < min_size){
      condition = F
      break
    }
    ## Set previous merge count so we can go back to it in case of premature exiting
    if(i == 1){
      prev_val = unaltered
    }else{prev_val = iterations[i - 1]}
    
    ## Make a merge vector to merge "val" adjacent points
    reduced_mappings$merge_vec <- floor(seq(from = 0, by = 1/val, length.out = nrow(reduced_mappings)))
    
    ## Merge data
    iteration_mappings <- reduced_mappings[,.(pos = first(pos), end = last(end), corrected_depth = sum(corrected_depth), n = length(pos), maf = merge_mafs(maf, exp = T), gc = mean(gc)), by = list(merge_vec, chr)]
    
    ## Map the previous iteration segment data onto current iteration data
    current_segment_mappings <- previous_segment_mappings[,-"end"][iteration_mappings, on = c("chr", "pos"), roll = Inf]
    current_segment_mappings[,segment_depth:=median(corrected_depth), by = list(chr, segment)]
    
    ## Previously flagged "copy-neutral" points used to compute the new "maxpeak" depth
    iteration_maxpeak <- median(maxpeak_segments[current_segment_mappings, on = c("chr", "segment")][(mp),]$corrected_depth)
    
    ## Readjust per-bin depth 
    current_segment_mappings[, corrected_depth := corrected_depth/(n) * val]
    current_segment_mappings[, n := NULL]
    current_segment_mappings[, segment_depth := median(corrected_depth), by = list(chr, segment)]
    
    iteration_positions <- subcl_pos/(unaltered/val)
    ## Adjust positions by the shift in maxpeak
    iteration_positions <- iteration_positions - (iteration_positions[names(iteration_positions) == ploidy] - iteration_maxpeak)
    ## Make an LM for regressing subclonal cnvs with the current set of depths
    df_depths = data.frame(cn = as.numeric(names(iteration_positions)), segment_depth = iteration_positions)
    lm_depths = lm(cn ~ segment_depth, data = df_depths)
    ## Get float CN values for per-chromosome adjustment of GMM means
    regress_cns <- round(predict(lm_depths, current_segment_mappings))
    ## Get per-chromosome CN values for per-chromosome adjustment of GMM means
    split_segs <- split(current_segment_mappings, f = current_segment_mappings$chr)
    individual_pos <- lapply(split_segs, function(x){
      vec <- as.numeric(names(table_vec(x$CN)))
      max_val <- max(ceiling(vec))
      ind_pos <- sort(c(vec, (0:max_val)[!0:max_val %in% vec]))
      cns <- names(table_vec(round(predict(lm_depths, x))))
      ind_pos <- sort(as.numeric(unique(c(ind_pos, cns))))
      ind_pos[ind_pos < 0] <- 0
      return(unique(ind_pos))
    })
    ## Get positions for this iteration
    iteration_positions <- depth(maxpeak = iteration_maxpeak, d = get_coverage_characteristics(tp, ploidy, iteration_maxpeak)$diff, P = ploidy, n = sort(unique(unlist(individual_pos))))
    ## Get a variance estimate of segments for this iteration
    iteration_var_nonmod <- weighted_median(current_segment_mappings[,.(sd_dp = sd(corrected_depth)), by = list(chr, segment)]$sd_dp, w = current_segment_mappings[,.(wt = length(corrected_depth)), by = list(chr, segment)]$wt)
    ## Get responsibilities based on current means/variance estimates
    wt = compute_responsibilities(current_segment_mappings$corrected_depth, iteration_positions, iteration_var_nonmod)
    ## Separate clonal and subclonal positions
    iteration_clonal_positions <- iteration_positions[which(as.numeric(names(iteration_positions)) == round(as.numeric(names(iteration_positions))))]
    iteration_subclonal_positions <- iteration_positions[which(as.numeric(names(iteration_positions)) != round(as.numeric(names(iteration_positions))))]
    ## Pick a variance estimate needed to separate the different clonal CNVs
    iteration_var = get_good_var(iteration_clonal_positions[1], iteration_clonal_positions[2], base = 5)
    ## get BAF variance per copy number
    mafstat <- current_segment_mappings[,.(vaf_sd = sd(unlist(unmerge_mafs(maf, flip = T))), n = length(na.omit(unmerge_mafs(maf)))), by = list(chr, segment)]
    maf_variance <- weighted.mean(mafstat$vaf_sd, mafstat$n, na.rm = T)
    ## Map CNs to segmented means based on new postion values
    current_segment_mappings$CN <- as.numeric(names(iteration_positions))[apply(parametric_gmm_fit(data = current_segment_mappings$segment_depth, means = iteration_positions, variances = iteration_var_nonmod), 1, which.max)]
    ## Get segment mappings
    collapsed_segs <- current_segment_mappings[,.(pos = first(pos), end = last(end), CN = first(CN), segment_depth = first(segment_depth), n = .N), by = list(chr, segment)]
    ## Another LM for CN from depth
    df_depths = data.frame(cn = as.numeric(names(iteration_positions)), segment_depth = iteration_positions)
    lm_depths = lm(cn ~ segment_depth, data = df_depths)
    ## Use current CN calls to get better fractional estimates per-chr
    collapsed_segs$fine_call = predict(lm_depths, collapsed_segs)
    current_segment_mappings$fine_call <- predict(lm_depths, current_segment_mappings)
    if(nrow(collapsed_segs) < nrow(subclonal_seg_mappings)){
      subclonal_seg_mappings <- collapsed_segs
    }
    collapsed_segs <- current_segment_mappings[,.(pos = first(pos), end = last(end), CN = first(CN), segment_depth = first(segment_depth), n = .N, scn = first(fine_call)), by = list(chr, segment)]
    pos_list <- lapply(1:nrow(collapsed_segs), function(x){
      cns = individual_pos[[collapsed_segs[x]$chr]]
      cns[cns == collapsed_segs[x]$CN] <- collapsed_segs[x]$scn
      cns
    })
    collapsed_segs <- collapsed_segs[,.(pos=first(pos), end = last(end), segment_depth = weighted.mean(segment_depth, w = n), n = sum(n)), by = list(chr, segment, CN)]
    collapsed_segs$segment <- 1:nrow(collapsed_segs)
    collapsed_segs$scn = predict(lm_depths, collapsed_segs)
    individual_pos <- lapply(unique(collapsed_segs$chr), function(x){
      chr_data = unique(round(collapsed_segs[chr == x]$scn))
      individual_pos[[x]] <- unique(sort(c(individual_pos[[x]], chr_data)))
      return(individual_pos[[x]])
    })
    names(individual_pos) <- unique(collapsed_segs$chr)
    pos_list <- lapply(1:nrow(collapsed_segs), function(x){
      cns = individual_pos[[collapsed_segs[x]$chr]]
      cns[cns == collapsed_segs[x]$CN] <- collapsed_segs[x]$scn
      cns
    })
    ## Now that the per-segment/chr noisiness is corrected, compute likelihoods for focal CNV detection
    ## We compute GMM likelihoods for segmented depth and raw depth, and look for discordant calls
    current_joint_probs <- maf_gmm_fit_subclonal_prior_segments(depth_data = current_segment_mappings$corrected_depth, vaf_data = current_segment_mappings$maf, chr_vec = current_segment_mappings$chr, means = iteration_positions, variances = iteration_var, ploidy = ploidy, maxpeak = iteration_maxpeak, maf_variances = maf_variance, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
    segged_joint_probs <- maf_gmm_fit_subclonal_prior_segments(depth_data = current_segment_mappings$segment_depth, vaf_data = current_segment_mappings$maf, chr_vec = current_segment_mappings$chr, means = iteration_positions, variances = iteration_var, ploidy = ploidy, maxpeak = iteration_maxpeak, maf_variances = maf_variance, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
    
    current_joint_probs <- current_joint_probs$jp_tbl
    current_joint_resps <- current_joint_probs/rowSums(current_joint_probs)
    
    segged_joint_probs <- segged_joint_probs$jp_tbl
    segged_joint_resps <- segged_joint_probs/rowSums(segged_joint_probs)
    
    ##Compress transitions by segment
    compressed_joint_resps <- data.table::copy(segged_joint_probs)
    compressed_joint_resps <- cbind(compressed_joint_resps, current_segment_mappings[,c("chr", "segment")])
    compressed_joint_resps <- compressed_joint_resps[, lapply(.SD, mean), by = list(chr, segment)]
    ## Mark ends of chromosomes
    chr_ends <- which(!na_or_true(shift(x = compressed_joint_resps$chr, type = "lead") == compressed_joint_resps$chr))
    ## Convert likelihoods into responsibilities, set NAs to 0 (don't match any GMM component)
    compressed_joint_resps <- compressed_joint_resps[,c("chr", "segment"):= NULL]
    compressed_joint_resps <- compressed_joint_resps/rowSums(compressed_joint_resps)
    compressed_joint_resps[which(is.na(compressed_joint_resps[,1])),] <- 0
    ## Get transition metrics for each meta-bin
    metric <- rowSums(abs(current_joint_resps - segged_joint_resps))
    metric[is.na(metric)] <- 2
    ## Try to get a threshold to introduce segment breaks
    shift_vec <- rowSums(abs(apply(compressed_joint_resps, 2, diff)))[-chr_ends]
    break_metric <- quantile(metric, prob = 0.99)
    
    ## Simulate the null hypothesis (constant-CN segment) to get a conservative threshold from this
    set.seed(42069)
    simulation_results = c()
    v_tbl = current_segment_mappings[,.(v = as.double(sd(corrected_depth)), s = as.double(sum(end - pos)), m = as.double(mean(corrected_depth))), by = c("CN")][order(CN)]
    ## Variance for the simulation
    sim_var = max(current_segment_mappings[,.(v = as.double(sd(corrected_depth)), s = as.double(sum(end - pos)), m = as.double(mean(corrected_depth))), by = c("CN")][order(CN)][s > quantile(s, 0.75)]$v)
    ## Variance is linearly proportional to copy number, so we use a linear model to get 
    ## variance estimates for each copy number
    if(nrow(v_tbl[s > quantile(s, 0.75)]) >= 3){
      lm_var = lm(v ~ CN, data = v_tbl[s > quantile(s, 0.75)], weights = s)
    }else{
      lm_var = lm(v ~ CN, data = data.table(CN = c(ploidy - 1, ploidy, ploidy + 1), v = rep(v_tbl[which.max(s)]$v, 3)))
    }
    ## Function to perform a round of MC simulation for a segment with no true breakpoints
    simulate_threshold = function(median_depth, segment_variance, model_variance, component_means, sim_size = 10000){
      simulated_seg = rnorm(sim_size, mean = median_depth, sd = segment_variance)
      simulated_segged = parametric_gmm_fit(rep(mean(simulated_seg), sim_size), means = component_means, variances = model_variance)
      simulated_segged = simulated_segged/rowSums(simulated_segged)
      simulated_nonsegged = parametric_gmm_fit(simulated_seg, means = component_means, variances = model_variance)
      zeroes = which(rowSums(simulated_nonsegged) == 0)
      simulated_nonsegged = simulated_nonsegged/rowSums(simulated_nonsegged)
      simulated_nonsegged[zeroes,] = 0
      simulated_thresh = rowSums(abs(simulated_segged - simulated_nonsegged))
      return(max(simulated_thresh[-which.max(simulated_thresh)]))
    }
    ## Get a separate transition threshold for each copy number
    thresh_by_cn = lapply(v_tbl$CN[v_tbl$CN == round(v_tbl$CN)], function(sim_cn){
      simulation_results = rep(NA, 10)
      for(z in 1:10){
        simulation_results[z] = simulate_threshold(iteration_clonal_positions[names(iteration_clonal_positions) == sim_cn], predict(lm_var, data.table(CN = sim_cn)), iteration_var, iteration_clonal_positions)
      }
      return(structure(median(simulation_results), names = sim_cn))
    })
    ## Table of thresholds by copy number
    thresh_by_cn = data.table("t" = unlist(thresh_by_cn), "CN" = as.numeric(names(unlist(thresh_by_cn))))
    ## Fill missing GMM components between the two tables
    if(ncol(segged_joint_resps) < ncol(current_joint_resps)){
      segged_joint_resps[,names(current_joint_resps)[!names(current_joint_resps) %in% names(segged_joint_resps)]] <- 0
    }
    ## Table of transition scores and CNs
    metric_df = data.table("metric" = metric, "CN" = round(current_segment_mappings$CN))
    ## Merge thresholds per CN with table of scores, identify putatitve segment breaks
    metric_df = thresh_by_cn[metric_df, on = "CN"]
    breaks = metric_df$metric >= metric_df$t
    ## Flag regions to break if they were identified above, or if they have a big shift in depth
    ## consistent with a 2-copy change (1-copy can introduce noise in very high resolution data)
    current_segment_mappings$flagged <- breaks | abs(current_segment_mappings$corrected_depth - current_segment_mappings$segment_depth) > get_coverage_characteristics(tp, ploidy, iteration_maxpeak)$diff*2
    ## use graph-based flagging/breaking to break & merge data back
    edge_vec <- c(1, rep(2:(nrow(current_segment_mappings)-1), each = 2), nrow(current_segment_mappings))
    g <- graph(edges = c(1, rep(2:(nrow(current_segment_mappings)-1), each = 2), nrow(current_segment_mappings)), directed = F)
    g <- set_vertex_attr(g, name = "chr", value = current_segment_mappings$chr)
    ## Does the chromosome match the lead chromosome? Use to find ends
    current_segment_mappings[, match := !na_or_true(chr == data.table::shift(chr, type = "lead"))]
    ## We should always break the ends of the segments
    chr_ends <- which(current_segment_mappings$match)
    segment_ends <- which(diff(current_segment_mappings$segment) != 0)
    segment_ends <- sort(c(segment_ends - 1, segment_ends, segment_ends + 1))
    if(any(segment_ends > max(E(g)))){
      segment_ends <- segment_ends[-which(segment_ends > max(E(g)))]
    }
    ## Get indices which have NA likelihoods so we always break these
    outlier_points <- which(current_segment_mappings$flagged)
    na_inds <- which(is.na(apply(current_joint_resps[outlier_points,], 1, sum)))
    na_inds <- unique(c(na_inds, which(is.na(apply(segged_joint_resps[outlier_points,], 1, sum)))))
    if(length(na_inds)){
      na_points <- outlier_points[na_inds]
      outlier_points <- outlier_points[-na_inds]
    }else{na_points <- c()}
    ## because subclones are a scourge, and allowing focal subclonal CNV detection introduces a ton of noise,
    ## we only allow new focal CNVs to be clonal
    alt_fits <- as.numeric(names(iteration_positions)[apply(current_joint_resps[outlier_points,], 1, which.max)])
    orig_fits <- as.numeric(names(iteration_positions)[apply(segged_joint_resps[outlier_points,], 1, which.max)])
    to_subclonal <- which(alt_fits != round(alt_fits) & (alt_fits != orig_fits))
    ## Flag meta-bins which changed a lot
    outlier_points <- c(outlier_points[-to_subclonal], na_points)
    outlier_edges <- unique(sort(c(outlier_points - 1, outlier_points)))
    outlier_edges <- outlier_edges[outlier_edges > 0 & outlier_edges < max(edge_vec) - 1]
    ## concatenate regions to break segments at
    to_delete <- unique(sort(c(chr_ends, segment_ends, outlier_edges)))
    ## Break segments
    g <- delete_edges(g, to_delete)
    ## Get broken segments
    busted_segs <- components(g)$membership
    busted_segment_mappings <- data.table::copy(current_segment_mappings)
    busted_segment_mappings$segment <- busted_segs
    busted_segment_mappings[, segment_depth := median(corrected_depth), by = segment]
    ## Get new estimates
    df_depths = data.frame(cn = as.numeric(names(iteration_positions)), segment_depth = iteration_positions)
    lm_depths = lm(cn ~ segment_depth, data = df_depths)
    ## Re-compute mixture model fit
    busted_jp_tbl <- maf_gmm_fit_subclonal_prior_segments(depth_data = busted_segment_mappings$segment_depth, vaf_data = busted_segment_mappings$maf, chr_vec = busted_segment_mappings$chr, means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
    ## Re-segment to enable broken ends to heal with old segments
    healed_segments <- ploidetect_prob_segmentator(prob_mat = busted_jp_tbl$jp_tbl, ploidy = ploidy, chr_vec = busted_segment_mappings$chr, seg_vec = busted_segment_mappings$segment, verbose = T, dist_vec = current_segment_mappings$segment_depth, subclones_discovered = T)
    busted_segment_mappings$segment <- healed_segments
    busted_segment_mappings[,segment_depth:=median(corrected_depth), by = list(chr, segment)]
    ## Remember all the checks for bad subclonal cnvs that we did earlier? They don't 100% work
    ## Bin shifts that were predicted to result in clonal CNVs end up resulting in subclonal CNVs
    ## which should not have been allowed to occur. So we need to fix that.
    ## Not sure if this should live here or not, or if a function is needed. I think in a previous 
    ## version it was. But I'll keep it as-is.
    repair_subcl_segs <- function(means, variances, maf_variances, maxpeak, tp, ploidy, cn_list, pos_list, seg_tbl, previous_jp_tbl, busted_segment_mappings){
      merge_dt <- busted_segment_mappings[,c("chr", "segment")]
      n_segs <- nrow(unique(merge_dt))
      healed_jp_tbl <- maf_gmm_fit_subclonal_prior_segments(depth_data = busted_segment_mappings$segment_depth, vaf_data = busted_segment_mappings$maf, chr_vec = busted_segment_mappings$chr, means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
      off_the_scale = which(rowSums(healed_jp_tbl$jp_tbl) == 0)
      full_orig_fits <- as.numeric(names(previous_jp_tbl)[apply(previous_jp_tbl, 1, which.max)])
      full_orig_fits <- unlist(cbind(full_orig_fits, merge_dt)[,.(full_orig_fits = median(full_orig_fits)), by = list(chr, segment)][,-c("chr", "segment")])
      full_new_fits <- as.numeric(names(healed_jp_tbl$jp_tbl)[apply(healed_jp_tbl$jp_tbl,1,which.max)])
      full_new_fits <- unlist(cbind(full_new_fits, merge_dt)[,.(full_new_fits = median(full_new_fits)), by = list(chr, segment)][,-c("chr", "segment")]
      ## Check for regions which transition from clonal to subclonal
      subcl_transition <- which(abs(full_orig_fits - full_new_fits) < 1 & abs(full_orig_fits - full_new_fits) != 0 & full_new_fits != round(full_new_fits))
      prec_clonal <- na_or_true(shift(full_new_fits)[subcl_transition] == round(shift(full_new_fits))[subcl_transition])
      # Succeeding
      suc_clonal <- na_or_true(shift(full_new_fits, type = "lead")[subcl_transition] == round(shift(full_new_fits, type = "lead"))[subcl_transition])
      ## Filter
      subcl_transition <- subcl_transition[prec_clonal & suc_clonal]
      ## Check that they are equal in CN
      subcl_transition <- subcl_transition[na_or_true(shift(full_new_fits)[subcl_transition] == shift(full_new_fits, type = "lead")[subcl_transition])]
      ### Now get the chr,seg combinations with bad segments
      bad_seg_maps <- unique(merge_dt)[subcl_transition]
      ## Merge segments using a graph
      g <- graph(edges = c(1, rep(2:(nrow(busted_segment_mappings)-1), each = 2), nrow(busted_segment_mappings)), directed = F)
      
      ## Basically repeat the earlier process of getting chr/segment ends
      segment_ends <- which(busted_segment_mappings$segment != shift(busted_segment_mappings$segment, type = "lead"))
      chr_ends <- which(busted_segment_mappings$chr != shift(busted_segment_mappings$chr, type = "lead"))
      tot_ends <- sort(unique(segment_ends, chr_ends))
      ends_tbl <- busted_segment_mappings[tot_ends]
      ends_tbl$ind <- 1:nrow(ends_tbl)
      ## Get new segment ends
      if(nrow(bad_seg_maps) > 0){
        filt_ends = ends_tbl[bad_seg_maps,, on = c("chr", "segment")]$ind
        filt_ends = filt_ends[!is.na(filt_ends)]
        new_ends <- tot_ends[-filt_ends]
      }else{
        new_ends <- tot_ends
      }
      new_ends <- sort(unique(c(new_ends, chr_ends)))
      ## Delete edges corresponding to segment ends, extract segment mappings
      g <- delete_edges(g, new_ends)
      busted_segment_mappings$segment <- components(g)$membership
      busted_segment_mappings[,segment_depth:=median(corrected_depth), by = list(chr, segment)]
      ## Get new GMM of likelihoods
      merge_jp_tbl <- maf_gmm_fit_subclonal_prior_segments(depth_data = busted_segment_mappings$segment_depth, vaf_data = busted_segment_mappings$maf, chr_vec = busted_segment_mappings$chr, means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
      busted_segment_mappings$call <- as.numeric(names(merge_jp_tbl$jp_tbl)[apply(merge_jp_tbl$jp_tbl, 1, which.max)])
      busted_segment_mappings$CN <- busted_segment_mappings$call
      ## Fix cases where same CN in consecutive segments, which shouldn't happen, but it does. 
      fixed_breaks <- which(busted_segment_mappings$chr != shift(busted_segment_mappings$chr, type = "lead") | busted_segment_mappings$CN != shift(busted_segment_mappings$CN, type = "lead"))
      g <- graph(edges = c(1, rep(2:(nrow(busted_segment_mappings)-1), each = 2), nrow(busted_segment_mappings)), directed = F)
      g <- delete_edges(g, fixed_breaks)
      ## Re-map segments, get new GMM of likelihoods
      busted_segment_mappings$segment <- components(g)$membership
      busted_segment_mappings[,segment_depth:=median(corrected_depth), by = list(chr, segment)]
      merge_jp_tbl <- maf_gmm_fit_subclonal_prior_segments(depth_data = busted_segment_mappings$segment_depth, vaf_data = busted_segment_mappings$maf, chr_vec = busted_segment_mappings$chr, means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs)
      ## Get newest calls with the filtered breakpoints
      busted_segment_mappings$call <- as.numeric(names(merge_jp_tbl$jp_tbl)[apply(merge_jp_tbl$jp_tbl, 1, which.max)])
      busted_segment_mappings[off_the_scale]$call <- round(predict(lm_depths, busted_segment_mappings[off_the_scale]))
      busted_segment_mappings$CN <- busted_segment_mappings$call
      f_n_segs <- nrow(unique(busted_segment_mappings[,c("chr", "segment")]))
      return(busted_segment_mappings)
    }
    ## Run the above function to repair errors made through subclone inclusion
    busted_segment_mappings <- repair_subcl_segs(means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos, pos_list = pos_list, seg_tbl = collapsed_segs, previous_jp_tbl = segged_joint_resps, busted_segment_mappings = busted_segment_mappings)
    ## Compute segment lengths to check if segment n50 has gone off a cliff due to oversegmentation
    seg_lens <- busted_segment_mappings[,.("pos" = first(pos), "end" = last(end)), by = list(chr, segment)][,.(diff=end-pos)]$diff
    seg_lens <- seg_lens[seg_lens > 0]
    previous_n50 <- current_n50
    previous_median_length <- current_median_length
    current_n50 <- n50_fun(seg_lens)
    current_median_length <- median(seg_lens)
    ## iterate
    i = i + 1
    ## 
    obs_pos <- as.numeric(names(table_vec(subcl_seg$CN)))
    obs_pos <- sort(c(obs_pos, (0:10)[!0:10 %in% obs_pos]))
    subcl_pos <- depth(maxpeak = maxpeak, d = get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)$diff, P = ploidy, n = obs_pos)
    
    ## Recompute maxpeak and the regions used to calculate it
    closeness <- abs(busted_segment_mappings$segment_depth - iteration_maxpeak)
    maxpeak_segments <- unique(busted_segment_mappings[which(closeness < diff(iteration_clonal_positions)[1]/2), c("chr", "segment")])
    maxpeak_segments$mp <- T
    ## Iteration housekeeping
    
    ## This is for when we reach maximum resolution & leave the loop "gracefully"
    if(i > length(iterations)){
      cn_positions = get_coverage_characteristics(tp, ploidy, iteration_maxpeak)$cn_by_depth
      out_seg_mappings <- data.table::copy(busted_segment_mappings)
      condition = F
      cn_positions = get_coverage_characteristics(tp, ploidy, iteration_maxpeak)$cn_by_depth
    }
    
    ## Exit case, where the contiguity drops too far. In this case we break the loop and output the previous mappings
    if(current_n50 < previous_n50/2){
      condition = F
      i = max(1, i - 2)
      val = iterations[i]
      out_seg_mappings <- data.table::copy(overseg_mappings)
    }else{overseg_mappings <- data.table::copy(busted_segment_mappings)}
    
    ## If all is good and we're continuing, save previous_segment_mappings to be the current (good)
    ## mappings
    if(condition){
      previous_segment_mappings <- data.table::copy(busted_segment_mappings)
      out_seg_mappings <- data.table::copy(previous_segment_mappings)
      previous_segment_mappings <- previous_segment_mappings[,.(pos = first(pos), CN = median(CN)), by = list(chr, segment)]
      cn_positions = get_coverage_characteristics(tp, ploidy, iteration_maxpeak)$cn_by_depth
    }
  }
  ## Call zygosity on output segments
  ## First get variance of BAFs
  out_maf_sd <- out_seg_mappings[,.(maf_var = sd(unmerge_mafs(maf, flip = T)), n = length(maf)), by = list(chr, segment)]
  maf_var <- weighted.mean(out_maf_sd$maf_var, w = out_maf_sd$n, na.rm = T)
<<<<<<< HEAD
  out_seg_mappings[,call:=as.numeric(call)]
  out_seg_mappings[,c("zygosity", "A", "B") := gmm_loh(maf, call, tp, ploidy, maf_var), by = list(chr, segment)]
=======
  ## Output numeric CNs
  out_seg_mappings[,call:=as.numeric(call)]
  ## Call LOH for all segments
  loh_calls <- out_seg_mappings[,.(zygosity = gmm_loh(maf, call, tp, ploidy, maf_var), call = first(call)), by = list(chr, segment)]
>>>>>>> Almost done R/probabilistic_segmentation.R
  
  ## Plotting states are from 0-8
  states <- c(0:8, 8)
  states = data.table(state = states)
  states$state_cn <- c(0:2, 2:3, 3:4, 4:5, 5)
  states$zygosity <- c(rep("HOM", times = 2), rep(c("HET", "HOM"), times = 4))
  
  out_seg_mappings$state_cn <- pmin(5, round(out_seg_mappings$call))
  
  out_seg_mappings <- states[out_seg_mappings, on = c("state_cn", "zygosity")]
  
  #loh_calls <- loh_calls[,(names(loh_calls) %in% c("chr", "segment", "state", "zygosity")), with = F]
  
  #setcolorder(loh_calls, c("chr", "segment", "state", "zygosity"))
  
  #out_seg_mappings <- loh_calls[out_seg_mappings, on = c("chr", "segment")]
  #out_seg_mappings[,c("state", "zygosity", "A", "B")]
  
  out_seg_mappings <- out_seg_mappings[,c("chr", "pos", "end", "segment", "corrected_depth", "segment_depth", "maf", "call", "state", "zygosity", "A", "B")]
  setnames(out_seg_mappings, "call", "CN")
  cytoband_path = Sys.glob("resources/*/cytobands.txt")[1]
  
  CN_calls <- split(out_seg_mappings, f = out_seg_mappings$chr)

  cna_plots <- list()
  
  for(i in 1:length(CN_calls)){
    cna_plots[i] <- list(plot_ploidetect(CN_calls[[i]], cn_positions, cytoband_path))
  }
  
  chrs = suppressWarnings(as.numeric(names(CN_calls)))
  sortedchrs = sort(chrs)
  chrs = c(sortedchrs, names(CN_calls)[is.na(chrs)])
  
  cna_plots = cna_plots[order(order(chrs))]
  
  CN_calls <- do.call(rbind.data.frame, CN_calls)
  
  segged_CN_calls <- CN_calls[,.(pos = first(pos), end = last(end), CN = first(CN), state = first(state), zygosity = first(zygosity), segment_depth = first(segment_depth), A = first(A), B = first(B)), by = list(chr, segment)]
  
  metadata = list(cn_positions = cn_positions)
  
  return(list("cna_plots" = cna_plots, "cna_data" = CN_calls, "segged_cna_data" = segged_CN_calls, "calling_metadata" = metadata))
}


gmm_loh <- function(in_mafs, CN, tp, ploidy, var){
  mafs <- unmerge_mafs(in_mafs, flip = T)
  if(length(CN) > 1){
    CN = CN[1]
  }
  if(CN <= 1.25){
    return(list("HOM", CN, 0))
  }
  if(length(mafs) == 0){
    A = round(CN)/2
    return(list("HET", A, CN - A))
  }
  pred_mafs <- testMAF_sc(CN, tp)
  pred_alleles <- as.numeric(names(pred_mafs))
  
  pred_mafs <- pred_mafs[!pred_alleles < round(floor(max(pred_alleles))/2)]
  fit <- parametric_gmm_fit(mafs, pred_mafs, var)
  resp <- fit/rowSums(fit)
  out_lik <- colSums(fit * resp)
  A = as.numeric(names(out_lik)[which.max(out_lik)])
  if(CN - as.numeric(names(out_lik)[which.max(out_lik)]) <= 0.25){
    return(list("HOM", A, CN - A))
  }else{
    return(list("HET", A, CN - A))
  }
}

one_segmenter <- function(in_list, in_models, bin_size = 100000){
  # First load model parameters from inputs
  maxpeak = in_list$maxpeak
  tp = in_models$tp[1]
  ploidy = in_models$ploidy[1]
  
  # Data
  binned_data = data.table(in_list$segmented_data)
  raw_data = data.table(in_list$all_data)
  
  # Get some parameters from those
  cov_char = get_coverage_characteristics(tp, ploidy, maxpeak)
  
  
  # Initial segmentation
  init_deviation = estimateVariance(to_seg = in_list$segmented_data$corrected_depth, size = 10, compress_iters = 10, folds = 100)
  
  transition_lik = zt_p(0, cov_char$diff/4, init_deviation)
  
  initial_segs = unlist(lapply(lapply(split(in_list$segmented_data$corrected_depth, f = in_list$segmented_data$chr), FUN = seed_compress, compress_iters = 1, transition_lik = transition_lik, var = init_deviation), function(x)x$segs))
  
  binned_data$segment <- initial_segs
  
  
  binned_data[,segment_depth:=median(corrected_depth), by = c("chr", "segment")]
  
  binned_data[,CN:=cn_from_dp(segment_depth, maxpeak, tp, ploidy)]
  
  
  
  
  binned_segs <- binned_data[,.(pos = first(pos), end = last(end), CN= first(CN), segment_depth = first(segment_depth), sd = sd(corrected_depth[2:(first(n)-1)]), n = first(n)), by = list(chr, segment)]
  
  
  chrt = 1
  plot_segments(binned_data[chr == chrt]$pos, binned_data[chr == chrt]$corrected_depth, binned_data[chr == chrt]$segment)
  
  segments = binned_data[chr == chrt]$segment
  to_seg = binned_data[chr == chrt]$corrected_depth
  
  segments <- paste0("contig_", segments)
  
  plot_segments(binned_data[chr == chrt]$pos, binned_data[chr == chrt]$corrected_depth, segment_repairer(to_seg, segments, maxpeak, cov_char))
  plot_segments(binned_data[chr == chrt]$pos, binned_data[chr == chrt]$corrected_depth, segment_breaker(to_seg, segments, maxpeak, cov_char))
  
  segments = segment_breaker(to_seg, segments, maxpeak, cov_char)
  segments <- paste0("contig_", segments)
  
  segments = segment_repairer(to_seg, segments, maxpeak, cov_char)
  
  plot_segments(binned_data[chr == chrt]$pos, binned_data[chr == chrt]$corrected_depth, segments)
  
  
  segment_repairer <- function(to_seg, segments, maxpeak, cov_char){
    ## Step 0: Estimate SD
    tbl = data.table("d" = to_seg, "s" = segments)
    tbl2 = data.table::copy(tbl)
    tbl = tbl[,.(d = mean(d), v = sd(d), n = .N), by = s]
    v = sqrt(weighted.mean(tbl$v^2, w = tbl$n, na.rm = T))
    tbl2[,v:=sd(d), by = s]
    v_lm = lm(v~d, data = tbl2)
    
    if(nrow(tbl) == 1){
      return(rep(1, tbl$n[1]))
    }
    ## Step 1: Merge segment ends
    seg_ends = grep("end", segments)
    if(length(seg_ends)){
      print("x")
    }
    
    ## Step 2: Merge adjacent segments if appropriate
    adjacent_gmm = function(tbl, v_lm){
      lapply(1:nrow(tbl), function(x){
        means = c(tbl$d[x - 1], tbl$d[x], tbl$d[x+1])
        if(x == 1){
          means = c(NA, means)
        }
        print(means)
        parametric_gmm_fit(tbl[x]$d, means, predict(v_lm, data.frame("d" = means)))
      }) %>% do.call(rbind, .)
    }
    
    train_val = adjacent_gmm(tbl = data.table("d" = c(maxpeak - cov_char$diff/4, maxpeak, maxpeak + cov_char$diff/4)), v_lm)[2,1]
    
    neighbor_p = apply(adjacent_gmm(tbl, v_lm)[,-2], 1, max, na.rm = T)
    
    which_neighbor = apply(adjacent_gmm(tbl, v_lm)[,-2], 1, which.max)
    
    neighbor_p <- which(neighbor_p >= train_val | tbl$n == 1)
    
    if(length(neighbor_p) == 0){
      return(rep(1:nrow(tbl), tbl$n))
    }
    
    merge_edges = unique(neighbor_p + which_neighbor[neighbor_p] - 2)
    
    lenseg = length(to_seg)
    g = graph(edges = c(1, rep(2:(lenseg - 1), each = 2), lenseg), directed = F)
    seg_ends = cumsum(tbl$n)[-merge_edges]
    seg_ends = seg_ends[-length(seg_ends)]
    segments = components(delete_edges(g, seg_ends))$membership
    return(segments)
  }
  
  segment_breaker <- function(to_seg, segments, maxpeak, cov_char){
    ## Step 0: Estimate SD
    tbl = data.table("d" = to_seg, "s" = segments)
    tbl2 = data.table::copy(tbl)
    tbl = tbl[,.(d = mean(d), v = sd(d), n = .N), by = s]
    v = sqrt(weighted.mean(tbl$v^2, w = tbl$n, na.rm = T))
    tbl2[,v:=sd(d), by = s]
    v_lm = lm(v~d, data = tbl2)
    ## Step 1: Identify shifts over CN = 1
    tbl[,c:=cn_from_dp(d, maxpeak, cov_char$tp, cov_char$ploidy)]
    
    
    
    transition_gmm = function(tbl, to_seg, v_lm){
      means = c()
      lapply(1:nrow(tbl), function(x){
        means = c(tbl$d, tbl$d[x], tbl$d[x+1])
        if(x == 1){
          means = c(NA, means)
        }
        parametric_gmm_fit(tbl[x]$d, means, predict(v_lm, data.frame("d" = means)))
      }) %>% do.call(rbind, .)
    }
    
    
    break_points = which(abs(cn_from_dp(to_seg, maxpeak, cov_char$tp, cov_char$ploidy) - rep(tbl$c, tbl$n)) > 1)
    
    if(length(break_points) == 0){
      return(segments)
    }
    
    ## Convert vertex IDs to edge IDs
    lenseg = length(to_seg)
    break_points <- c(break_points, break_points - 1)
    break_points <- break_points[break_points > 0 & break_points < lenseg]
    
    break_points = unique(sort(c(break_points, cumsum(tbl$n))))
    break_points = break_points[-length(break_points)]
    
    ## Calc new segments
    g = graph(edges = c(1, rep(2:(lenseg - 1), each = 2), lenseg), directed = F)
    return(components(delete_edges(g, break_points))$membership)
    
  }
  
  bin_size = 100000
  reduced_data = data.table::copy(binned_data)
  
  {
    binned_segs = reduced_data[,.(pos = first(pos), end = last(end)), by = c("chr", "segment")]
    bin_size = bin_size/2
    
    reduced_data = ploidetect_preprocess(raw_data, simplify = T, simplify_size = bin_size)
    
    red_maxpeak = reduced_data$maxpeak
    cov_char_r = get_coverage_characteristics(tp, ploidy, red_maxpeak)
    
    reduced_data = reduced_data$x
    
    reduced_data = binned_segs[,c("chr", "pos", "segment")][reduced_data[,c("chr", "pos", "end", "corrected_depth", "maf")], on = c("chr", "pos"), roll = Inf]
    
    segments = unlist(lapply(split(reduced_data, f = reduced_data$chr), function(x){segment_breaker(to_seg = x$corrected_depth, segments = x$segment, red_maxpeak, cov_char_r)}))
    reduced_data$segment = segments
    segments = unlist(lapply(split(reduced_data, f = reduced_data$chr), function(x){segment_repairer(to_seg = x$corrected_depth, segments = x$segment, red_maxpeak, cov_char_r)}))
    
    reduced_data$segment = segments
    
    chr_t = 1
    plot_segments(reduced_data[chr == chr_t]$pos, reduced_data[chr == chr_t]$corrected_depth, reduced_data[chr == chr_t]$segment)
    
    
  }
  
  reduced_data[,CN:=cn_from_dp(corrected_depth, red_maxpeak, tp, ploidy)]
  reduced_data %>% filter(chr == 1) %>% ggplot(aes(x = pos, y = corrected_depth, color = CN)) + geom_point() + scale_color_viridis()
  
  
  # Initial segmentation
  init_deviation = estimateVariance(to_seg = in_list$segmented_data$corrected_depth, size = 10, compress_iters = 10, folds = 100)
  
  transition_lik = zt_p(0, cov_char$diff/4, init_deviation)
  
  
  
  binned_data[,n:=.N, by = c("chr", "segment")]
  
  binned_segs <- binned_data[,.(pos = first(pos), end = last(end), CN= first(CN), segment_depth = first(segment_depth), sd = sd(corrected_depth[2:(first(n)-1)]), n = first(n)), by = list(chr, segment)]
  
  binned_segs[,parent_cn:=round(CN)]
  
  binned_segs[,size:=end-pos]
  
  subcl_test_segs <- binned_segs[CN < ploidy + 3]
  
  subcl_test_segs[,fraction:=abs(CN - parent_cn)]
  
  subcl_test_segs[order(fraction, decreasing = T)]
  
  subcl_test_segs <- split(subcl_test_segs, subcl_test_segs$parent_cn)
  
  for(i in 1:nrow(subcl_test_segs)){
    x = unlist(subcl_test_segs[1])
    ks.test(binned_data[chr == x[1] & segment == x[2]], rnorm(as.numeric(x[8]), mean = depth(maxpeak, cov_char$diff, ploidy, as.numeric(x[9]))))$p.value
  }
  
  apply(subcl_test_segs, 1, function(x){
    print(x)
    ks.test(binned_data[chr == x[1] & segment == x[2]]$corrected_depth, rnorm(10000, mean = depth(maxpeak, cov_char$diff, ploidy, as.numeric(x[9])), sd = as.numeric(x[7])))$p.value
  })
}

plot_segments <- function(chr_vec, chr, pos, y, segments){
  pos = pos[chr_vec == chr]
  y = y[chr_vec == chr]
  segments = segments[chr_vec == chr]
  seg_pos <- pos[segments != shift(segments)]
  p = data.frame(pos, y) 
  p = p %>% ggplot(aes(x= pos, y = y, color = segments)) + geom_point() + geom_vline(xintercept = seg_pos, alpha=0.1) + scale_color_viridis()
  return(p)
}
