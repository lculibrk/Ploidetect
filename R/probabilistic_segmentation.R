ploidetect_prob_segmentator <- function(prob_mat, ploidy, chr_vec, seg_vec, dist_vec, verbose = T, lik_shift = 1, subclones_discovered = F){
  if(verbose){
    print("Performing segmentation and copy number variant calling")
  }
  if(!all(c("chr", "seg") %in% names(prob_mat))){
    prob_mat <- cbind(prob_mat, "chr" = chr_vec, "seg" = seg_vec)
  }
  
  orig_mat <- data.table::copy(prob_mat)
  
  resp_mat <- data.table::copy(prob_mat)
  
  resp_mat <- resp_mat[,-((ncol(resp_mat)-3):ncol(resp_mat))]
  
  resp_mat <- resp_mat/rowSums(resp_mat)
  
  resp_mat[which(is.na(resp_mat[,1])),] <- 0
  
  compressed_mat <- resp_mat[, lapply(.SD, mean), by = list(orig_mat$chr, orig_mat$seg)]
  
  compressed_mat <- compressed_mat[,-c(1:2)]
  if(nrow(compressed_mat) < nrow(orig_mat)){
    states <- apply(compressed_mat, 1, which.max)
    state_transitions <- data.frame("leading" = states[-1], "lagging" = states[1:(length(states)-1)])
    state_transitions <- split(state_transitions, f = 1:nrow(state_transitions))
    state_transitions <- lapply(state_transitions, unlist)
    
    shifted <- lagged_df(compressed_mat)
    shifted <- shifted[-nrow(shifted),]
    
    #transition_liks <- lapply(1:length(state_transitions), function(x){
    #  sum(abs(compressed_mat[x,state_transitions[[x]], with = F]/sum(compressed_mat[x,state_transitions[[x]], with = F]) - shifted[x, state_transitions[[x]], with = F]/sum(shifted[x, state_transitions[[x]], with = F])))
    #})
    
    #transition_liks <- unlist(transition_liks)
    #transition_liks <- transition_liks[!na_or_true(transition_liks == 0)]
    
    all_transitions <- apply(compressed_mat, 2, diff) %>% abs %>% apply(1, sum)
    transition_threshold <- min(1.5, quantile(all_transitions[all_transitions > 0.5], prob = 0.5))
  }else{transition_threshold = lik_shift}
  
  if(is.null(nrow(prob_mat))){
    print(prob_mat)
  }
  prob_mat$chr <- chr_vec
  prob_mat$seg <- seg_vec
  prob_mat$dist <- dist_vec
  prob_mat[which(is.na(prob_mat[,1])),] <- 0
  datsplit <- split(as.data.frame(prob_mat), prob_mat$chr)
  ## Run segmentation by compression on all chromosomes
  if(verbose){
    print("Performing segmentation of copy number data")
  }
  compressedalldat <- unlist(lapply(datsplit, runiterativecompression_prob, segmentation_threshold = transition_threshold, verbose = verbose, subclones_discovered = subclones_discovered))
  return(compressedalldat)
}

runiterativecompression_prob <- function(data, segmentation_threshold = segmentation_threshold, verbose = F, subclones_discovered = F){

  if(verbose){
    print("Running iterative compression to segment read depth data")
  }
  segments <- data$seg
  dist = data$dist
  
  data <- data[,-which(names(data) %in% c("chr", "seg", "dist"))]
  message(paste0("Initial segment count: ", length(unique(segments))))
  converged <- F
  compress <- data.table(data)
  compress <- split(compress, f = segments)
  if(length(compress) == 1){
    converged = T
  }


  while(!converged){
    windows <- length(compress)
    
    compress <- compressdata_prob(compress, criteria = segmentation_threshold, dist_vec = dist, subclones_discovered = subclones_discovered)


    if(length(compress) == windows | length(compress) == 1){
      converged <- T
    }
  }
  segs <- rep(1:length(compress), times = lapply(compress, nrow))
  return(segs)
}

compressdata_prob <- function(compress, criteria, dist_vec, subclones_discovered = F){
  reps <- lapply(compress, nrow)
  ids = rep(x = 1:length(reps), times = unlist(reps))
  
  dists_original <- dist_vec
  
  dists <- aggregate(dist_vec, list(ids), FUN = mean)
  
  compress_original = rbindlist(compress)
  
  compressed_compress <- compress_original[, lapply(.SD, mean), by = factor(ids)][,-1]
  
  rel_liks <- data.table::copy(compressed_compress)
  
  #rel_liks <- compressed_compress/rowSums(compressed_compress)
  
  #rel_liks[apply(rel_liks, 1, function(x)all(is.na(x))),] <- 0
  # Arrange data by position and remove tibble-ness
  # Get differences between neighbouring points

  if(subclones_discovered){
    states <- apply(rel_liks, 1, which.max)
    state_probs <- apply(rel_liks, 1, max)
    state_transitions <- data.frame("leading" = states[-1], "lagging" = states[1:(length(states)-1)])
    state_transitions <- split(state_transitions, f = 1:nrow(state_transitions))
    state_transitions <- lapply(state_transitions, unlist)
    shifted <- lagged_df(rel_liks)
    t_liks <- as.data.frame(rel_liks)
    
    state_vec <- 1:length(states)
    
    
    
    
    fit_val <- vapply(1:length(states), function(x){t_liks[x,states[x]]}, 0.01)
    fit_next <- vapply(1:length(states), function(x){t_liks[x,c(states, 1)[x+1]]}, 0.01)
    
    fit_df <- cbind(fit_val, fit_next)
    fit_df <- fit_df/rowSums(fit_df)
    
    fit_shifted <- vapply(1:length(states), function(x){t_liks[x + 1, states[x]]}, 0.01)
    fit_shifted_next <- vapply(1:length(states), function(x){t_liks[x+1,c(states, 1)[x+1]]}, 0.01)
    
    fit_shifted_df <- cbind(fit_shifted, fit_shifted_next)
    fit_shifted_df <- fit_shifted_df/rowSums(fit_shifted_df)
    
    
    
    transition_liks <- rowSums(abs(fit_df - fit_shifted_df))
    
    transition_probs <- transition_liks[-length(transition_liks)]
    transition_probs[which(is.na(transition_probs))] <- dists$x[which(is.na(transition_probs))]
  }else{
    transition_probs <- rowSums(abs(apply(rel_liks, 2, diff)))
  }
  #transition_prob <- sapply(1:nrow(state_transitions), function(x){
  #  transition_probs[x, state_transitions$leading[x]] - transition_probs[x, state_transitions$lagging[x]]
  #})

  
  #shifted - rel_liks
  
  #shift(rel_liks)
  
  #transition_probs <- apply(rel_liks, 2, diff)
  #if(!is.null(nrow(transition_probs))){
  #  transition_probs <- rowSums(abs(transition_probs))
  #}else{transition_probs <- sum(abs(transition_probs))}

  
  
  #if(length(compress) > 2){
  #  transition_probs <- apply(transition_probs, 1, function(x)sum(abs(x)))
  #}else if(length(compress) == 2){
  #  transition_probs <- max(transition_probs)
  #}
  # Initialize graph from data
  graph <- graph(edges = c(row.names(compressed_compress)[1], rep(row.names(setDF(compressed_compress)[-c(1, nrow(compressed_compress)),]), each = 2), row.names(setDF(compressed_compress))[nrow(setDF(compressed_compress))]), directed = F)
  # Set edges to have the diffs attribute
  graph <- set_edge_attr(graph, name = "transition", value = transition_probs)
  # Give vertices appropriate attributes
  graph <- set_vertex_attr(graph, name = "probs", value = compress)
  graph <- set_vertex_attr(graph, name = "npoints", value = unlist(reps))
  # loop over all vertices
  sort_by_diff <- data.frame("vertex" = V(graph)$name)
  sort_by_diff$diff <- 0
  #for(row in 1:nrow(sort_by_diff)){
  #  sort_by_diff$diff[row] <- max(edge_attr(graph, "diff", incident(graph, sort_by_diff$vertex[row])))
  #  sort_by_diff$e[row] <- which.max(edge_attr(graph, "diff", incident(graph, sort_by_diff$vertex[row])))
  #}
  edges <- edge_attr(graph, "transition", E(graph))
  e1 <- c(NA, edges)
  e2 <- c(edges, NA)
  
  diff_vec <- abs(diff(dist_vec))
  d1 <- shift(dist_vec, type = c("lag"))
  d2 <- shift(dist_vec, type = c("lead"))
  
  #new_seg_data %>% filter(chr == "X") %>%  mutate(segment = ids) %>%  filter() %>% ggplot(aes(x = pos, y = corrected_depth, color = states)) + geom_point() + scale_color_viridis()
  
  e_comp <- as.numeric(na_or_true(e1 > e2))
  e_eq <- which(e1 == e2)
  d_comp <- as.numeric(na_or_true(d1 > d2))
  
  e_comp[e_eq] <- d_comp[e_eq]
  
  if(all(unlist(reps) == 1)){
    preserve <- unique(1:length(V(graph)) - (1-e_comp))
    sequentials <- preserve[which(na_or_true((preserve - shift(preserve) == 1)))]
    sequentials <- sequentials[!sequentials %in% c(1, length(edges) + 1)]
    preserve <- preserve[!preserve %in% (sequentials - e_comp[sequentials])]
    to_del <- V(graph)[!V(graph) %in% preserve]
    graph <- delete_edges(graph, to_del)
  }else{
    ## Diff check
    # Get incident vertices of cases where both edges are beyond the threshold
    both_ind <- na_or_true(e1 > criteria) & na_or_true(e2 > criteria)
    both_ind <- which(both_ind & (unlist(reps) > 1))
    both <- unique(sort(pmax(1, c(both_ind, both_ind - 1))))
    
    both <- both[both < length(e1)]
    
    ## Check for vertices where there is only one bin in the segment
    lenone <- which(unlist(reps) == 1)
    ## Preserve the edge that has the lower differential probability
    preserve <- unique(lenone - (1 - e_comp[lenone]))
    
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
    
    to_delete <- unique(c(both, single))
    to_delete <- to_delete[!to_delete %in% preserve]
    
    
    graph <- delete_edges(graph, to_delete)
  }
  #na.omit(unique((1:length(e1)) - as.numeric(e1 > e2)))
  

  #sort_by_diff <- sort_by_diff %>% arrange(diff)
  #time26 <- Sys.time()
  #delete_edges(graph, E(graph)[edge_attr(graph, "diff") > segmentation_threshold*x])
  #todel <- c()
  #for(vertex in sort_by_diff$vertex){
  #  if(length(incident(graph, vertex)) == 0){
  #    next
  #  }
  #  # If vertex is an outlier (diffs are over some threshold fraction of what we expect for a copy change) then break all edges
  #  if(all(edge_attr(graph, "diff", incident(graph, vertex)) > segmentation_threshold*x)){
  #    todel <- c(todel, incident(graph, vertex))
  #    #graph <- delete_edges(graph, incident(graph, vertex))
  #    next
  #  }
  #  # If vertex has two edges, break the one with larger "diff"
  #  #if(length(incident(graph, vertex)) == 2){
  #  #  graph <- delete_edges(graph, incident(graph, vertex)[which.max(get.edge.attribute(graph, "diff", incident(graph, vertex)))])
  #  #}
  #  if(length(incident(graph, vertex)) == 2){
  #    todel <- c(todel, incident(graph, vertex)[which.max(get.edge.attribute(graph, "diff", incident(graph, vertex)))])
  #  }
  #}
  #time27 <- Sys.time()
  #todel <- E(graph)[edge_attr(graph, "diff") > segmentation_threshold*x]
  #graph <- delete_edges(graph, todel)
  #delete_edges(graph, todel)
  #time3 <- Sys.time()
  #time3-time26
  #time27-time26
  #print(get.vertex.attribute(toy, "npoints", V(toy)))
  #print(toy)
  
  # Get list of all vertex pairs to merge
  tomerge <- ends(graph, E(graph))
  merging <- as.numeric(c(t(tomerge)))
  diff(merging)[1:(length(merging) - 1) %% 2 == 0]
  
  
  
  #groups(graph)
  # Get all vertices
  #vertices <- V(graph)
  # Coerce vertices into a format where value is the vertex value and name is vertex name
  #vertnames <- names(vertices)
  #vertices <- as.numeric(vertices)
  #names(vertices) <- vertnames
  # Change "tomerge" from names to values
  #tomerge[,2] <- vertices[which(names(vertices) %in% tomerge[,2])]
  #tomerge[,1] <- vertices[which(names(vertices) %in% tomerge[,1])]
  # Not needed I think
  #todelete <- vertices[which(vertices %in% tomerge[,2])]
  # Change pairs of vertices to repeat the same vertex twice (used in contract.vertices() to map which vertices to contract)
  #vertices[which(vertices %in% tomerge[,2])] <- tomerge[,1]
  #mode(vertices) <- "integer"
  #lint <- vertices[1]
  #time4 <- Sys.time()
  #for(i in 2:length(vertices)){
  #  if(vertices[i] - 1 > lint){
  #    vertices[vertices == vertices[i]] <- lint + 1
  #  }
  #  lint <- vertices[i]
  #}
  #time5 <- Sys.time()
  
  # Merge connected vertices
  
  vertices <- components(graph)$membership
  
  graph <- contract.vertices(graph, mapping = vertices, vertex.attr.comb = list("probs" = function(...){
    elements = c(...)
    #elements = lapply(elements, unlist)
    return(rbindlist(elements))
  }, "npoints" = "sum", "name" = "first"))
  
  
  # Delete all old vertices
  #toy <- delete.vertices(toy, which(names(V(toy)) == "character(0)"))
  # Reconstruct the data.frame we began with
  out_probs = get.vertex.attribute(graph, "probs")
  n = get.vertex.attribute(graph, "npoints")
  
  compress <- out_probs
  i_segs <- 1:length(compress)
  i_segs <- rep(i_segs, times = unlist(lapply(compress, nrow)))
  subcl_seg %>% filter(chr == "X") %>% mutate("seg" = i_segs, "breakpoint" = i_segs != shift(i_segs)) %>%  filter(pos < 0.5e+07) %>% ggplot(aes(x = pos, y = corrected_depth, color = seg == 6)) + geom_point() + scale_color_viridis(discrete = T)
  subcl_seg %>% filter(chr == "X") %>% mutate("seg" = i_segs, "breakpoint" = i_segs != shift(i_segs)) %>%  filter() %>% ggplot(aes(x = pos, y = corrected_depth, color = seg)) + geom_point() + scale_color_viridis(discrete = F)
  
  #print(dat[which(dat$npoints == 1),])
  if(T){
    message(paste0("Iteration complete, segment count = ", length(compress)))
  }
  #current_segment_mappings %>% filter(chr == 11) %>%  mutate(segment = i_segs) %>%  filter(pos < 1e+08, pos > 5e+07) %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis()
  
  #print(c(time2 - time1, time2-time3, time4-time3, time5-time4, time6-time5))
  #print(dat)
  return(compress)
}

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
  
  ## Correct depth
  unmerged_data$corrected_depth <- unmerged_data$corrected_depth + unmerged_maxpeak
  
  em_sd <- match_kde_height(data = unmerged_data$corrected_depth, means = predictedpositions, sd = em_sd)
  
  props <- parametric_gmm_fit(unmerged_data$corrected_depth, means = predictedpositions, variances = em_sd)
  props <- colSums(props)/sum(colSums(props))
  
  pdf_fun <- mixpdf_function(predictedpositions, props, sd = em_sd)
  
  
  
  den <- density(unmerged_data$corrected_depth[unmerged_data$corrected_depth < max(predictedpositions)], n = 2^16)
  
  den_df <- data.frame("x" = den$x, "dens" = den$y, "pred" = pdf_fun(den$x)$y)
  
  den_df %>% ggplot(aes(x = x, y = dens)) + geom_line() + geom_line(data = den_df, aes(x = x, y = pred, color = "pdf"))
  
  
  maf_gmm_result <- maf_gmm_fit(depth_data = unmerged_data$corrected_depth, vaf_data = unmerged_data$maf, chr_vec = unmerged_data$chr, means = predictedpositions, variances = em_sd, maf_variances = segment_maf_sd, ploidy = ploidy, maxpeak = unmerged_maxpeak)
  ## segment
  
  
  
  #segs <- ploidetect_prob_segmentator(prob_mat = as.matrix(maf_gmm_result$jp_tbl), ploidy = 2, chr_vec = unmerged_data$chr)
  
  
  
  
  #unmerged_data$segment <- segs
  
  train_df <- data.frame("CN" = as.numeric(names(means)), "segment_depth" = means)
  train <- lm(CN ~ segment_depth, train_df)
  data.table(unmerged_data)[,.(segment_depth = median(corrected_depth)), by = list(chr, segment)]
  
  unmerged_data <- unmerged_data %>% group_by(chr, segment) %>% dplyr::mutate("segment_depth" = median(corrected_depth))
  
  unmerged_data$CN = round(predict(train, unmerged_data))
  
  #unmerged_data %>% filter(chr == 13) %>% ggplot(aes(x = pos, y = corrected_depth, color = CN)) + geom_point() + scale_color_viridis()

  
  calls <- cbind(unmerged_data$chr, unmerged_data$segment, maf_gmm_result$jp_tbl)
  calls <- calls[, lapply(.SD, median), by = list(V1, V2)]
  calls$CN <- names(calls)[-(1:2)][apply(calls, 1, function(x)which.max(x[-(1:2)]))]

  names(calls)[1:2] <- c("chr", "segment")
  
  calls <- calls[,c("chr", "segment", "CN")]
  

  unmerged_data <-   left_join(unmerged_data, calls, by = c("chr", "segment"))
  #unmerged_data %>% filter(chr == 12) %>% ggplot(aes(x = pos, y = corrected_depth, color = as.numeric(apply(compressed_compress, 1, max)))) + geom_point() + theme(legend.position = "None") + scale_color_viridis(discrete = F)
  
  
  
  calls <- aggregate(data = as.matrix(maf_gmm_result$jp_tbl), FUN = "median")
  calls <- data.frame("CN" = names(calls)[calls[,-1] %>% apply(1, which.max)], segment = calls$Group.1)
  
  
  ## Sanitise highoutliers and merge with rest of data
  #unmerged_highoutliers <- unmerged_highoutliers[,c("tumour", "normal", "maf", "wind", "size", "gc", "tumour", "size")]
  #names(unmerged_highoutliers) <- names(unmerged_data$x)

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
  sd <- grouped_data %>% group_by(chr, segment) %>% dplyr::summarise("sd" = sd(corrected_depth), "mean_residual" = median(corrected_depth), "length" = n()) %>% ungroup %>% arrange(desc(length)) %>% slice(n()/2) %>% summarise("medsd" = median(sd, na.rm = T)) %>% unlist
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
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% dplyr::mutate("mean_residual" = median(corrected_depth), "z" = (corrected_depth - median(corrected_depth))/sd)
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

subclonal_joint_probs <- function(new_seg_data, tp, ploidy, maxpeak, clonal_positions, vaf_variance, fracs = c(0.25, 0.5), in_sd = NA){
  predictedpositions <- get_coverage_characteristics(tp, ploidy, maxpeak)$cn_by_depth
  ## Generate an LM for regressing subclonal CNAs by depth
  train_df <- data.frame("segment_depth" = clonal_positions, "CN" = as.numeric(names(predictedpositions)))
  model <- lm(CN ~ segment_depth, train_df)
  ## Get depth differential from TP/ploidy estimate
  d_diff <- get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)$diff
  
  ## Pick which fraction we will use for subclonal CNA calling
  int_liks <- c()
  variances <- c()
  positions_vectors <- list()
  for(subcl_frac in fracs){
    ## Call subclonal CNAs by position as an estimate
    new_seg_data$CN <- round(predict(model, new_seg_data)/(10*subcl_frac), digits = 1)*10*subcl_frac
    
    ## Limit subclonal CNA calling to ploidy+3 tops
    subcl_pos <- seq(from = min(round(new_seg_data$CN)), to = max(round(new_seg_data$CN)), by = subcl_frac)
    ## Generate a vector of depths
    subcl_pos <- structure(depth(maxpeak, d = d_diff, P = ploidy, n = seq(from = min(round(new_seg_data$CN)), to = max(round(new_seg_data$CN)), by = subcl_frac)), names = subcl_pos)
    to_add <- (0:10)[which(!0:10 %in% names(subcl_pos))]
    
    subcl_pos <- sort(c(subcl_pos, depth(maxpeak, d = d_diff, P = ploidy, n = to_add)))
    subcl_pos <- subcl_pos[-which(as.numeric(names(subcl_pos)) > ploidy + 3 & as.numeric(names(subcl_pos)) - floor(as.numeric(names(subcl_pos))) != 0)]
    
    ## Store positions
    positions_vectors <- c(positions_vectors, list(subcl_pos))
    if(is.na(in_sd)){
      ## Compute variance by segment, grab only largest 50% of segments or top ten, whichever is smaller
      segged_sd <- new_seg_data %>% 
        filter(CN - round(CN) == 0) %>%  
        group_by(chr, segment) %>% 
        add_tally %>% 
        dplyr::summarise("var" = sd(corrected_depth), "n" = first(n)) %>% 
        arrange(desc(n)) %>% 
        ungroup %>%  
        slice(1:min(10, nrow(.)/2)) %>% 
        select(var) %>% 
        unlist %>% 
        mean
      ## Use the previous estimate to get the true variance by comparison to a KDE
      segged_sd <- match_kde_height(new_seg_data$segment_depth, means = subcl_pos, sd = segged_sd)
    }else{segged_sd <- in_sd}
    ## Store variance
    variances <- c(variances, segged_sd)
    
    ## Fit segmented depth to new subclonal GMM
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
  subcl_pos <- positions_vectors[[which.max(int_liks)]]
  segged_sd <- variances[[which.max(int_liks)]]
  print(segged_sd)
  

  
  #subclonal_probs <- maf_gmm_fit_subclonal(depth_data = new_seg_data$segment_depth, chr_vec = new_seg_data$chr, vaf_data = new_seg_data$maf, means = subcl_pos, variances = segged_sd, maf_variances = vaf_variance, maxpeak = maxpeak, ploidy = ploidy)
  subclonal_probs <- data.table(parametric_gmm_fit(data = new_seg_data$segment_depth, means = subcl_pos, variances = segged_sd))
  colnames(subclonal_probs) <- names(subcl_pos)
  subclonal_probs <- list("jp_tbl" = subclonal_probs)
  #plot_density_gmm(data = new_seg_data$segment_depth, means = subcl_pos[1:30], weights = colSums(subclonal_probs$jp_tbl, na.rm = T)/sum(colSums(subclonal_probs$jp_tbl, na.rm = T)), sd = segged_sd)
  
  return(list("fraction" = fracs[which.max(int_liks)], "probs" = subclonal_probs, "segged_sd" = segged_sd))
}

segment_subclones <- function(new_seg_data, predictedpositions, depth_variance, vaf_variance, maxpeak, tp, ploidy){
  
  ## Generate an LM for regressing subclonal CNAs by depth
  train_df <- data.frame("segment_depth" = predictedpositions, "CN" = as.numeric(names(predictedpositions)))
  model <- lm(CN ~ segment_depth, train_df)
  
  ## Compute joint probabilities
  subclonal_probs <- subclonal_joint_probs(new_seg_data = new_seg_data, tp = tp, ploidy = ploidy, maxpeak = maxpeak, clonal_positions = predictedpositions, vaf_variance = 0.03)
  
  ##Store fraction
  subclonal_fraction <- subclonal_probs$fraction
  subclonal_variance <- subclonal_probs$segged_sd
  subclonal_probs <- subclonal_probs$probs

  
  ## Debugging plot
  
  plot_density_gmm(data = new_seg_data$segment_depth, means = depth(maxpeak, d = get_coverage_characteristics(tp, ploidy, maxpeak)$diff, P = ploidy, n = as.numeric(names(subclonal_probs$jp_tbl))), weights = colSums(subclonal_probs$jp_tbl)/sum(colSums(subclonal_probs$jp_tbl)), sd = subclonal_variance)
  
  ## Generate segments from subclonal probability matrix using signal compression segmentation
  newvelle_segs <- ploidetect_prob_segmentator(prob_mat = subclonal_probs$jp_tbl, ploidy = ploidy, chr_vec = new_seg_data$chr, seg_vec = unlist(lapply(split(1:nrow(new_seg_data), new_seg_data$chr), function(x)1:length(x))), dist_vec = new_seg_data$corrected_depth)
  
  ## Map segments
  #new_seg_data <- new_seg_data %>% mutate("segment" = newvelle_segs, "prob" = probs) %>% group_by(chr, segment) %>% dplyr::mutate(segment_depth = mean(corrected_depth), prob = median(prob))
  
  new_seg_data$segment <- newvelle_segs
  
  new_seg_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  
  #new_seg_data %>% filter(chr == 5, pos < 3e+07, pos > 2e+07, segment_depth < 2e+05) %>%  ggplot(aes(x = pos, y = corrected_depth, color = segment == 14)) + geom_point() + scale_color_viridis(discrete = T)
  
  new_seg_data$CN <- round(predict(model, data.frame("segment_depth" = new_seg_data$segment_depth))/(subclonal_fraction*10), 1)*(subclonal_fraction*10)
  
  CNA_list <- new_seg_data[,c("chr", "CN")]
  
  CNA_list <- unique(CNA_list)

  
  CNA_list$CN[CNA_list$CN > ploidy + 3] <- round(CNA_list$CN[CNA_list$CN > ploidy + 3])
  comp_pos <- depth(maxpeak = maxpeak, d = get_coverage_characteristics(tp, ploidy, maxpeak)$diff, P = ploidy, sort(unique(CNA_list$CN)))
  
  CNA_list$CN <- pmax(0, CNA_list$CN)
  
  CNA_list <- split(CNA_list$CN, f = CNA_list$chr)
  CNA_list <- lapply(CNA_list, function(x)unique(sort(x)))
  
  #new_jp_tbl <- maf_gmm_fit_subclonal_prior(depth_data = new_seg_data$segment_depth, vaf_data = new_seg_data$maf, chr_vec = new_seg_data$chr, means = comp_pos, variances = depth_variance, maf_variances = 0.06, maxpeak = maxpeak, ploidy = ploidy, tp = tp, cn_list = CNA_list)
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
  
  ## Fix subclonal segments that might be noisy
  edge_vec <- c(1, rep(2:(nrow(new_seg_data)-1), each = 2), nrow(new_seg_data))
  g <- graph(edges = edge_vec, directed = F)
  
  g <- set_vertex_attr(g, name = "chr", value = new_seg_data$chr)
  
  todel <- c(which(abs(diff(new_seg_data$parent_cns)) > 0), 
             which(!na_or_true(shift(new_seg_data$chr, type = "shift") == new_seg_data$chr)))
  g <- delete_edges(g, todel)
  
  new_seg_data$segment <- components(g)$membership
  ## Debugging plot
  new_seg_data %>% mutate(fit = names(subclonal_probs$jp_tbl)[apply(subclonal_probs$jp_tbl, 1, which.max)]) %>%  filter(chr == "9") %>% ggplot(aes(x = pos, y = segment_depth, color = parent_cns)) + geom_point() + scale_color_viridis()
  
  return(list("data" = list(new_seg_data), "fraction" = subclonal_fraction, "subclonal_variance" = subclonal_variance))
}

#' @export
ploidetect_cna_sc <- function(all_data, segmented_data, tp, ploidy, maxpeak, verbose = T, min_size = 1){
  predictedpositions <- get_coverage_characteristics(tp, ploidy, maxpeak)$cn_by_depth
  
  
  variance <- density(segmented_data$corrected_depth)$bw
  variance <- match_kde_height(segmented_data$corrected_depth, means = predictedpositions, sd = variance)

  proportions <- compute_responsibilities(segmented_data$corrected_depth, means = predictedpositions, variances = variance)
  proportions <- colSums(proportions)/sum(colSums(proportions))
  
  ploidy <- as.numeric(names(proportions)[which.max(proportions)])
  maxpeak <- predictedpositions[which.max(proportions)]
  
  variance <- match_kde_height(segmented_data$corrected_depth, means = predictedpositions, sd = variance, comparison_point = maxpeak)
  
  
  
  #joint_probs <- maf_gmm_fit(depth_data = segmented_data$corrected_depth, vaf_data = segmented_data$maf, chr_vec = segmented_data$chr, means = predictedpositions, variances = variance, ploidy = ploidy, maxpeak = maxpeak, maf_variances = 0.06)
  
  joint_probs <- list("jp_tbl" = data.table(parametric_gmm_fit(segmented_data$corrected_depth, means = predictedpositions, variances = variance)))
  plot_density_gmm(data = segmented_data$corrected_depth, means = depth(maxpeak, diff(predictedpositions)[1], ploidy, as.numeric(names(joint_probs$jp_tbl))), weights = colSums(joint_probs$jp_tbl/rowSums(joint_probs$jp_tbl, na.rm = T), na.rm = T), sd = variance)
  clonal_cnas <- ploidetect_prob_segmentator(prob_mat = joint_probs$jp_tbl, ploidy = ploidy, chr_vec = segmented_data$chr, seg_vec = unlist(lapply(split(1:nrow(segmented_data), segmented_data$chr), function(x)1:length(x))), dist_vec = segmented_data$corrected_depth, lik_shift = 0.75)
  clonal_cna_data <- segmented_data
  clonal_cna_data$segment <- clonal_cnas
  gc_tbl <- data.table(all_data[,c("chr", "pos", "gc")])
  gc_tbl$size <- all_data$end - all_data$pos
  clonal_cna_data <- data.table(clonal_cna_data)
  clonal_cna_data$row = 1:nrow(clonal_cna_data)
  gc_tbl <- clonal_cna_data[,c("chr", "pos", "row")][gc_tbl, on = list(chr, pos), roll = Inf]
  gc_tbl <- gc_tbl[,.(gc = mean(gc), size = sum(size)), by = row]
  gc_tbl <- gc_tbl[!is.na(row)]
  clonal_cna_data$gc <- gc_tbl$gc
  clonal_cna_data$size <- gc_tbl$size
  
  

  clonal_cna_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  fit <- maf_gmm_fit(depth_data = segmented_data$segment_depth, vaf_data = segmented_data$maf, chr_vec = segmented_data$chr, means = predictedpositions, variances = variance, ploidy = ploidy, maxpeak = maxpeak, maf_variances = 0.06)
  clonal_cna_data$fit <- apply(joint_probs$jp_tbl, 1, which.max)
  
  clonal_cna_data[,corrected_depth:= segment_depth + lowesswrapper(gc, corrected_depth, bw = 0.75)$residual, by = list(chr, segment)]
  clonal_cna_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  
  joint_probs <- list("jp_tbl" = data.table(parametric_gmm_fit(clonal_cna_data$corrected_depth, means = predictedpositions, variances = variance)))
  clonal_cnas <- ploidetect_prob_segmentator(prob_mat = joint_probs$jp_tbl, ploidy = ploidy, chr_vec = segmented_data$chr, seg_vec = unlist(lapply(split(1:nrow(segmented_data), segmented_data$chr), function(x)1:length(x))), dist_vec = segmented_data$corrected_depth, lik_shift = 0.75)
  clonal_cna_data$segment <- clonal_cnas
  clonal_cna_data[,segment_depth := median(corrected_depth), by = list(chr, segment)]
  
  clonal_cna_data %>% filter(chr == "X") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis(discrete = F) + geom_point(aes(x = pos, y = segment_depth))
  #clonal_cna_data %>% filter(fit == 10) %>% ggplot(aes(x = gc, y =corrected_depth, color = segment)) + geom_point() + scale_color_viridis(discrete = F)# + geom_line(aes(x = pos, y = gc * maxpeak*2)) + scale_y_continuous(sec.axis = sec_axis(trans = ~./(maxpeak*2)))
  
  #clonal_cna_data %>% filter(chr == "6") %>% ggplot(aes(x = gc, y = corrected_depth)) + geom_point() + geom_smooth(method = "loess", span = 1)
  
  subcl_seg <- segment_subclones(new_seg_data = clonal_cna_data, predictedpositions = predictedpositions, depth_variance = variance, vaf_variance = 0.06, maxpeak = maxpeak, tp = tp, ploidy = ploidy)
  subclonal_fraction <- subcl_seg$fraction
  subclonal_variance <- subcl_seg$subclonal_variance
  subcl_seg <- subcl_seg$data[[1]]
  subcl_seg$subclonal <- (as.numeric(subcl_seg$CN) - round(subcl_seg$CN)) != 0
  
  common_call <- names(which.max(table_vec(subcl_seg$CN)))
  common_maf_means <- testMAF(as.numeric(common_call), tp)
  maf_variance <- match_kde_height(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, sd = 0.03)
  
  t <- colSums(compute_responsibilities(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, variances = maf_variance))
  
  plot_density_gmm(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, sd = maf_variance, weights = t)
  
  ## Subclonal positions
  obs_pos <- as.numeric(names(table_vec(subcl_seg$CN)))
  obs_pos <- sort(c(obs_pos, (0:10)[!0:10 %in% obs_pos]))

  subcl_seg$CN[subcl_seg$CN < 0] <- 0
  
  subcl_cn <- as.numeric(names(table_vec(subcl_seg$CN)))
  subcl_pos <- depth(maxpeak, d_diff, ploidy, subcl_seg$CN)
  
  subcl_seg$dev_pos <- subcl_seg$corrected_depth - subcl_pos
  
  subcl_seg %>% filter(chr == "X") %>% ggplot(aes(x = pos, y = corrected_depth, color = corrected_depth - n_corrected_depth)) + geom_point() + scale_color_viridis(discrete = F)
  subcl_seg %>% filter(chr == "X") %>% ggplot(aes(x = gc, y = dev_pos)) + geom_point() + geom_smooth(method = "loess", span = 0.5)
  
  subcl_seg[,corrected_depth := lowesswrapper(x = gc, y = dev_pos, bw = 0.5)$residual + segment_depth, by = list(chr, segment)]
  
  clonal_cna_data$corrected_depth <- subcl_seg$corrected_depth
  
  subcl_seg <- segment_subclones(new_seg_data = clonal_cna_data, predictedpositions = predictedpositions, depth_variance = variance, vaf_variance = 0.06, maxpeak = maxpeak, tp = tp, ploidy = ploidy)
  subclonal_fraction <- subcl_seg$fraction
  subclonal_variance <- subcl_seg$subclonal_variance
  subcl_seg <- subcl_seg$data[[1]]
  subcl_seg$subclonal <- (as.numeric(subcl_seg$CN) - round(subcl_seg$CN)) != 0
  
  common_call <- names(which.max(table_vec(subcl_seg$CN)))
  common_maf_means <- testMAF(as.numeric(common_call), tp)
  maf_variance <- match_kde_height(as.numeric(unlist(unmerge_mafs(subcl_seg$maf[subcl_seg$CN == common_call]))), means = common_maf_means, sd = 0.03)
  
  #subcl_seg %>% filter(chr == "X") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis(discrete = F)
  
  
  
  split_segs <- split(subcl_seg, f = subcl_seg$chr)
  
  model <- lm(CN ~ corrected_depth, data = data.frame(CN = as.numeric(names(predictedpositions)), corrected_depth = predictedpositions))
  
  individual_pos <- lapply(split_segs, function(x){
    vec <- as.numeric(names(table_vec(x$CN)))
    max_val <- max(ceiling(vec))
    ind_pos <- sort(c(vec, (0:max_val)[!0:max_val %in% vec]))
    cns <- names(table_vec(round(predict(model, x))))
    ind_pos <- sort(as.numeric(unique(c(ind_pos, cns))))
    ind_pos[ind_pos < 0] <- 0
    return(unique(ind_pos))
  })
  
  unaltered <- ploidetect_preprocess(all_data, simplify = T, simplify_size = 100000, verbose = F)$merged
  
  all_data_preprocessed <- ploidetect_preprocess(all_data, simplify = F, simplify_size = NA, verbose = T)
  all_data_preprocessed <- all_data_preprocessed$x
  
  initial_segment_mappings <- subcl_seg[,.(pos = first(pos), CN = first(CN)), by = list(chr, segment)]
  
  reduced_mappings <- initial_segment_mappings[all_data_preprocessed, on = c("chr", "pos"), roll = Inf]
  
  reduced_mappings[, segment_depth := median(corrected_depth), by = list(chr, segment)]
  
  #d_red <- density(reduced_mappings$segment_depth, n = 2^16)
  
  #reduced_maxpeak <- d_red$x[which.max(d_red$y)]
  
  #p_dp <- depth(reduced_maxpeak, d = get_coverage_characteristics(tp, ploidy, reduced_maxpeak)$diff, P = ploidy, n = reduced_mappings$CN)
  
  #reduced_mappings$dev_pos <- reduced_mappings$corrected_depth - p_dp
  
  #reduced_mappings %>% filter(chr == "X") %>% ggplot(aes(x = pos, y = dev_pos, color = segment)) + geom_point() + scale_color_viridis()
  #reduced_mappings %>% filter(chr == "X", segment == 263) %>% ggplot(aes(x = gc, y = dev_pos, color = segment)) + geom_point() + scale_color_viridis() + geom_smooth(method = "loess", span = 0.75)
  
  #approx_bins <- lapply(unique(subcl_seg$chr), function(x)seq(from = min(subcl_seg[chr == x]$pos), to = max(subcl_seg[chr == x]$end), by = 100000))
  #names(approx_bins) <- unique(subcl_seg$chr)
  
  #approx_bins <- data.table(stack(approx_bins))
  #names(approx_bins) <- c("pos", "chr")
  #approx_bins$meta_bin <- 1:nrow(approx_bins)
  
  #reduced_mappings <- approx_bins[reduced_mappings, on = c("chr", "pos"), roll = Inf]
  
  #head(reduced_mappings)
  
  #reduced_mappings[, t_corrected_depth := lowesswrapper(gc, dev_pos, 1)$residual + segment_depth, by = list(meta_bin, segment)]
  #gcfit = function(x, y){
  #  if(length(x) >=3){
  #    return(lm(y ~ poly(x, 3))$residual)
  #  }
  #  else{
  #    return(lm(y ~ x)$residual)
  #  }
  #}
  
  reduced_mappings %>% filter(chr == "4") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis() #+ scale_y_continuous(limits = c(0, 5000))
  reduced_mappings %>% filter(chr == "4", segment == 158) %>% ggplot(aes(x = gc, y = corrected_depth)) + geom_point(size = 5, alpha = 0.5) + geom_smooth(method = "loess", span = 800)
  
  
  #plot_density_gmm(data = reduced_mappings$segment_depth, means = base_characteristics$cn_by_depth, weights = rep(1, times = 11), sd = variance/unaltered)
  
  iterations <- round(unaltered/2^(1:ceiling(log2(unaltered))))
  
  ## Re-estimate maxpeak
  ploidy <- as.numeric(names(subcl_pos)[which.min(abs(density(subcl_seg$segment_depth, n = 2^16)$x[which.max(density(subcl_seg$segment_depth, n = 2^16)$y)] - subcl_pos))])
  maxpeak <- density(subcl_seg$segment_depth, n = 2^16)$x[which.max(density(subcl_seg$segment_depth, n = 2^16)$y)]
  
  closeness <- abs(subcl_seg$segment_depth - maxpeak)
  maxpeak_segments <- unique(subcl_seg[which(closeness < diff(predictedpositions)[1]/2), c("chr", "segment")])
  
  maxpeak_segments$mp <- T
  
  #subcl_seg %>% ggplot(aes(x = pos, y = segment_depth, color = closeness < diff(predictedpositions)[1]/2)) + geom_point()
  
  subcl_pos <- depth(maxpeak = maxpeak, d = get_coverage_characteristics(tp = tp, ploidy = ploidy, maxpeak = maxpeak)$diff, P = ploidy, n = obs_pos)
  
  maxpeak_base <- maxpeak/(unaltered - 1)
  
  base_characteristics <- get_coverage_characteristics(tp, ploidy, maxpeak_base)
  
  previous_segment_mappings <- data.table::copy(initial_segment_mappings)
  
  seg_lens <- diff(previous_segment_mappings[,.("pos" = first(pos)), by = list(chr, segment)]$pos)
  seg_lens <- seg_lens[seg_lens > 0]
  current_n50 <- n50_fun(seg_lens)
  current_median_length <- median(seg_lens)
  print("iter")
  condition = T
  i = 1
  while(condition){
    
    val <- iterations[i]
    
    if(val < min_size){
      condition = T
      break
    }
    
    if(i == 1){
      prev_val = unaltered
    }else{prev_val = iterations[i - 1]}
    
    reduced_mappings$merge_vec <- floor(seq(from = 0, by = 1/val, length.out = nrow(reduced_mappings)))
    
    iteration_mappings <- reduced_mappings[,.(pos = first(pos), end = last(end), corrected_depth = sum(corrected_depth), n = length(pos), maf = merge_mafs(maf, exp = T), gc = mean(gc)), by = list(merge_vec, chr)]
    
    current_segment_mappings <- previous_segment_mappings[iteration_mappings, on = c("chr", "pos"), roll = Inf]
    
    current_segment_mappings %>% filter(chr == "4") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis()
    
    iteration_maxpeak <- median(maxpeak_segments[current_segment_mappings, on = c("chr", "segment")][(mp),]$corrected_depth)
    
    current_segment_mappings[, corrected_depth := corrected_depth/(n) * val]
    current_segment_mappings[, n := NULL]
    current_segment_mappings[, segment_depth := median(corrected_depth), by = list(chr, segment)]
    
    current_segment_mappings[,corrected_depth:= lowesswrapper(x = gc, y = corrected_depth, bw = 1)$residual + median(corrected_depth), by = list(CN)]
    
    #current_segment_mappings %>% filter(chr == "17") %>% ggplot(aes(x = pos, y = t_cor_dp, color = segment)) + geom_point() + scale_color_viridis()
    #current_segment_mappings %>% filter(chr == "17") %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis()
    
    current_segment_mappings %>% group_by(chr, segment) %>% dplyr::summarise(v = sd(corrected_depth)) %>% ungroup %>% select(v) %>% unlist %>% mean(na.rm = T)
    current_segment_mappings %>% group_by(chr, segment) %>% dplyr::summarise(v = sd(t_cor_dp)) %>% ungroup %>% select(v) %>% unlist %>% mean(na.rm = T)
    
    
    #current_segment_mappings %>% filter(chr == "13", segment == 56) %>% ggplot(aes(x = gc, y = t_cor_dp, color = segment)) + geom_point() + scale_color_viridis() + geom_smooth(method = "loess", span = 1)
    
    
    iteration_positions <- subcl_pos/(unaltered/val)
    
    iteration_positions <- iteration_positions - (iteration_positions[names(iteration_positions) == ploidy] - iteration_maxpeak)
    
    
    
    iteration_var <- weighted_median(current_segment_mappings[,.(sd_dp = sd(corrected_depth)), by = list(chr, segment)]$sd_dp, w = current_segment_mappings[,.(wt = length(corrected_depth)), by = list(chr, segment)]$wt)
    #iteration_var <- match_kde_height(current_segment_mappings$corrected_depth, iteration_positions, iteration_var)
    
    
    
    wt = compute_responsibilities(current_segment_mappings$corrected_depth, iteration_positions, iteration_var)
    
    plot_density_gmm(data = current_segment_mappings$corrected_depth, means = iteration_positions, weights = colSums(wt), sd = iteration_var)
    
    
    iteration_clonal_positions <- iteration_positions[which(as.numeric(names(iteration_positions)) == round(as.numeric(names(iteration_positions))))]
    iteration_subclonal_positions <- iteration_positions[which(as.numeric(names(iteration_positions)) != round(as.numeric(names(iteration_positions))))]
    
        

    
    current_segment_mappings %>% filter(chr == "11", corrected_depth, corrected_depth < max(iteration_positions)) %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis(discrete = F)
    
    
    mafstat <- current_segment_mappings[,.(vaf_sd = sd(unlist(unmerge_mafs(maf, flip = T))), n = length(na.omit(unmerge_mafs(maf)))), by = list(chr, segment)]
    

    
    maf_variance <- weighted.mean(mafstat$vaf_sd, mafstat$n, na.rm = T)
    
    current_joint_probs <- maf_gmm_fit_subclonal_prior(depth_data = current_segment_mappings$corrected_depth, vaf_data = current_segment_mappings$maf, chr_vec = current_segment_mappings$chr, means = iteration_positions, variances = iteration_var, ploidy = ploidy, maxpeak = iteration_maxpeak, maf_variances = maf_variance, tp = tp, cn_list = individual_pos)
    #fill_cols <- data.table(matrix(0, ncol = length(iteration_subclonal_positions), nrow = nrow(current_joint_probs$jp_tbl)))
    #colnames(fill_cols) <- names(iteration_subclonal_positions)
    segged_joint_probs <- maf_gmm_fit_subclonal_prior(depth_data = current_segment_mappings$segment_depth, vaf_data = current_segment_mappings$maf, chr_vec = current_segment_mappings$chr, means = iteration_positions, variances = iteration_var, ploidy = ploidy, maxpeak = iteration_maxpeak, maf_variances = maf_variance, tp = tp, cn_list = individual_pos)

    
    
    current_joint_probs <- current_joint_probs$jp_tbl
    current_joint_resps <- current_joint_probs/rowSums(current_joint_probs)
    
    segged_joint_probs <- segged_joint_probs$jp_tbl
    segged_joint_resps <- segged_joint_probs/rowSums(segged_joint_probs)
    
    compressed_joint_resps <- data.table::copy(segged_joint_probs)
    compressed_joint_resps <- cbind(compressed_joint_resps, current_segment_mappings[,c("chr", "segment")])
    compressed_joint_resps <- compressed_joint_resps[, lapply(.SD, mean), by = list(chr, segment)]
    chr_ends <- which(!na_or_true(shift(x = compressed_joint_resps$chr, type = "lead") == compressed_joint_resps$chr))
    compressed_joint_resps <- compressed_joint_resps[,c("chr", "segment"):= NULL]
    #seg_metric <- quantile(apply(compressed_joint_resps, 1, max), probs = 0.25)
    compressed_joint_resps <- compressed_joint_resps/rowSums(compressed_joint_resps)
    compressed_joint_resps[which(is.na(compressed_joint_resps[,1])),] <- 0
    
    metric <- rowSums(abs(current_joint_resps - segged_joint_resps))
    metric[is.na(metric)] <- 2
    
    #metric[which(apply(current_joint_probs, 1, max) <= seg_metric)] <- 0
    
    break_metric <- quantile(rowSums(abs(apply(compressed_joint_resps, 2, diff)))[-chr_ends], prob = 0.5)
    
    current_segment_mappings %>% mutate("prob" = rowSums(current_joint_probs)) %>%  filter(corrected_depth < max(iteration_positions), chr == "6") %>% ggplot(aes(x = pos, y = corrected_depth, color = prob)) + geom_point() + scale_color_viridis(discrete = F) + geom_hline(yintercept = iteration_positions[1:20])
    
    
    
    if(ncol(segged_joint_resps) < ncol(current_joint_resps)){
      segged_joint_resps[,names(current_joint_resps)[!names(current_joint_resps) %in% names(segged_joint_resps)]] <- 0
    }

    
    current_segment_mappings %>% mutate("prob" = metric) %>%  filter(chr == "1") %>% ggplot(aes(x = pos, y = corrected_depth, color = prob)) + geom_point() + scale_color_viridis(discrete = F)
    current_segment_mappings %>% mutate("prob" = metric) %>%  filter(chr == "13") %>% ggplot(aes(x = pos, y = corrected_depth, color = apply(segged_joint_resps, 1, which.max))) + geom_point() + scale_color_viridis(discrete = F)
    current_segment_mappings %>% mutate("prob" = metric) %>%  filter(chr == "13") %>% ggplot(aes(x = pos, y = corrected_depth, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    
    
    current_segment_mappings$flagged <- metric >= break_metric
    edge_vec <- c(1, rep(2:(nrow(current_segment_mappings)-1), each = 2), nrow(current_segment_mappings))
    g <- graph(edges = c(1, rep(2:(nrow(current_segment_mappings)-1), each = 2), nrow(current_segment_mappings)), directed = F)
    g <- set_vertex_attr(g, name = "chr", value = current_segment_mappings$chr)
    
    current_segment_mappings[, match := !na_or_true(chr == data.table::shift(chr, type = "lead"))]
    
    chr_ends <- which(current_segment_mappings$match)
    segment_ends <- which(diff(current_segment_mappings$segment) != 0)
    segment_ends <- sort(c(segment_ends - 1, segment_ends, segment_ends + 1))
    if(any(segment_ends > max(E(g)))){
      segment_ends <- segment_ends[-which(segment_ends > max(E(g)))]
    }
    

    outlier_points <- which(current_segment_mappings$flagged)
    na_inds <- which(is.na(apply(current_joint_resps[outlier_points,], 1, sum)))
    na_inds <- unique(c(na_inds, which(is.na(apply(segged_joint_resps[outlier_points,], 1, sum)))))
    if(length(na_inds)){
      na_points <- outlier_points[na_inds]
      outlier_points <- outlier_points[-na_inds]
    }else{na_points <- c()}

    
    alt_fits <- as.numeric(names(iteration_positions)[apply(current_joint_resps[outlier_points,], 1, which.max)])
    orig_fits <- as.numeric(names(iteration_positions)[apply(segged_joint_resps[outlier_points,], 1, which.max)])
    to_subclonal <- which(alt_fits != round(alt_fits) & (alt_fits != orig_fits))
    
    outlier_points <- c(outlier_points[-to_subclonal], na_points)
    outlier_edges <- unique(sort(c(outlier_points - 1, outlier_points)))
    outlier_edges <- outlier_edges[outlier_edges > 0 & outlier_edges < max(edge_vec) - 1]
    
    to_delete <- unique(sort(c(chr_ends, segment_ends, outlier_edges)))
    
    print("del_edges")
    
    print(unique(sort(chr_ends)))
    print(unique(sort(segment_ends)))
    print(unique(sort(outlier_edges)))
    
    g <- delete_edges(g, to_delete)
    
    busted_segs <- components(g)$membership
    busted_segment_mappings <- data.table::copy(current_segment_mappings)
    busted_segment_mappings$segment <- busted_segs
    busted_segment_mappings[, segment_depth := median(corrected_depth), by = segment]
    
    
    
    #busted_segment_mappings %>% mutate("prob" = apply(current_joint_probs, 1, max) < seg_metric) %>%  filter(chr == "X", segment_depth < max(iteration_positions)) %>% ggplot(aes(x = pos, y = corrected_depth, color = prob)) + geom_point() + scale_color_viridis(discrete = T) + geom_hline(yintercept = iteration_clonal_positions[iteration_clonal_positions < 15000])
    
    busted_jp_tbl <- maf_gmm_fit_subclonal_prior(depth_data = busted_segment_mappings$segment_depth, vaf_data = busted_segment_mappings$maf, chr_vec = busted_segment_mappings$chr, means = iteration_positions, variances = variance/(unaltered/val), maf_variances = maf_variance, maxpeak = iteration_maxpeak, ploidy = ploidy, tp = tp, cn_list = individual_pos)
    
    healed_segments <- ploidetect_prob_segmentator(prob_mat = busted_jp_tbl$jp_tbl, ploidy = ploidy, chr_vec = busted_segment_mappings$chr, seg_vec = busted_segment_mappings$segment, verbose = T, dist_vec = current_segment_mappings$segment_depth)
    
    busted_segment_mappings$segment <- healed_segments
    
    busted_segment_mappings %>% mutate("prob" = metric) %>%  filter(chr == "6", segment_depth < max(iteration_positions)) %>% ggplot(aes(x = pos, y = corrected_depth, color = segment)) + geom_point() + scale_color_viridis(discrete = F)
    seg_lens <- diff(busted_segment_mappings[,.("pos" = first(pos)), by = list(chr, segment)]$pos)
    seg_lens <- seg_lens[seg_lens > 0]
    previous_n50 <- current_n50
    previous_median_length <- current_median_length
    current_n50 <- n50_fun(seg_lens)
    current_median_length <- median(seg_lens)
    cn_calls <- data.table::copy(busted_jp_tbl$jp_tbl)
    cn_calls <- cbind("chr" = busted_segment_mappings$chr, "segment" = busted_segment_mappings$segment, cn_calls)
    cn_calls <- cn_calls[,lapply(.SD, mean), .SDcols = 3:ncol(cn_calls), by = list(chr, segment)]
    cn_calls[,call:=as.numeric(names(cn_calls[,3:ncol(cn_calls)])[which.max(.SD)]), by = 1:nrow(cn_calls), .SDcols = 3:ncol(cn_calls)]
    cn_calls <- cn_calls[,c(1:2, ncol(cn_calls)), with = F]
    busted_segment_mappings <- cn_calls[busted_segment_mappings, on = c("chr", "segment")]
    busted_segment_mappings %>% mutate("prob" = apply(current_joint_probs, 1, max)) %>% filter(chr == "6", segment_depth < max(iteration_positions)) %>% ggplot(aes(x = pos, y = corrected_depth, color = call)) + geom_point() + scale_color_viridis(discrete = F)
    
    i = i + 1
    if(i > length(iterations)){
      condition = F
    }
    if(current_n50 < previous_n50/2){
      condition = F
      i = max(1, i - 2)
      val = iterations[i]
      
    }
    if(condition){
      previous_segment_mappings <- data.table::copy(busted_segment_mappings)
      out_seg_mappings <- data.table::copy(previous_segment_mappings)
      previous_segment_mappings <- previous_segment_mappings[,.(pos = first(pos)), by = list(chr, segment)]
    }
  }
  
  out_seg_mappings %>% filter(chr == "1", segment_depth < max(subcl_pos/(unaltered/val)), call < 6) %>% ggplot(aes(x = pos, y = corrected_depth, color = call)) + geom_point() + scale_color_viridis(discrete = F)
  
  out_maf_sd <- out_seg_mappings[,.(maf_var = sd(unmerge_mafs(maf, flip = T)), n = length(maf)), by = list(chr, segment)]
  maf_var <- weighted.mean(out_maf_sd$maf_var, w = out_maf_sd$n, na.rm = T)
  loh_calls <- out_seg_mappings[,.(zygosity = gmm_loh(maf, call, tp, ploidy, maf_var), call = first(call)), by = list(chr, segment)]
  
  states <- c(0:8, 8)
  states = data.table(state = states)
  states$state_cn <- c(0:2, 2:3, 3:4, 4:5, 5)
  states$zygosity <- c(rep("HOM", times = 2), rep(c("HET", "HOM"), times = 4))
  
  loh_calls$state_cn <- pmin(5, round(loh_calls$call))
  
  loh_calls <- states[loh_calls, on = c("state_cn", "zygosity")]
  
  loh_calls <- loh_calls[,(names(loh_calls) %in% c("chr", "segment", "state", "zygosity")), with = F]
  
  setcolorder(loh_calls, c("chr", "segment", "state", "zygosity"))
  
  
  out_seg_mappings <- loh_calls[out_seg_mappings, on = c("chr", "segment")]
  
  CN_palette <- c("0" = "#cc0000", 
                  "1" = "#000066", 
                  "2" = "#26d953", 
                  "3" = "#609f70", 
                  "4" = "#cccc00", 
                  "5" = "#80804d",
                  "6" = "#cc6600", 
                  "7" = "#856647", 
                  "8" = "#cc0000"
  )
  CN_calls <- split(out_seg_mappings, f = out_seg_mappings$chr)
  
  CNA_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)[1]] | pos > centromeres$end[which(centromeres$chr %in% chr)[2]]) %>% ggplot(aes(x = pos, y = log(corrected_depth + maxpeak), color = as.character(state))) + 
      geom_point(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("log(Read Depth)") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
      theme_bw()
  })
  vaf_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)][1] | pos > centromeres$end[which(centromeres$chr %in% chr)][2]) %>% filter(!is.na(maf)) %>% ggplot(aes(x = pos, y = unlist(unmerge_mafs_grouped(maf, flip = T)), color = as.character(state))) + 
      geom_point(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("Major allele frequency") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " allele frequency profile")) + 
      scale_y_continuous(limits = c(0.5, 1)) +
      theme_bw()
  })
  
  cna_plots <- list()
  
  for(i in 1:length(CNA_plot)){
    cna_plots[i] <- list(plot_grid(CNA_plot[[i]], vaf_plot[[i]], align = "v", axis = "l", ncol = 1))
  }
  
  CN_calls <- do.call(rbind.data.frame, CN_calls)
  
  
  return(list("cna_plots" = cna_plots, "cna_data" = CN_calls))
}


gmm_loh <- function(in_mafs, CN, tp, ploidy, var){
  mafs <- unmerge_mafs(in_mafs, flip = T)
  if(length(CN) > 1){
    CN = CN[1]
  }
  if(CN < 2){
    return("HOM")
  }
  if(length(mafs) == 0){
    return("HET")
  }
  pred_mafs <- testMAF_sc(CN, tp)
  pred_alleles <- as.numeric(names(pred_mafs))
  pred_mafs <- pred_mafs[-which(pred_alleles < round(floor(max(pred_alleles))/2))]
  fit <- parametric_gmm_fit(mafs, pred_mafs, var)
  resp <- fit/rowSums(fit)
  out_lik <- colSums(fit * resp)
  if(floor(as.numeric(names(out_lik)[which.max(out_lik)])) == floor(CN)){
    return("HOM")
  }else{
    return("HET")
  }
}

