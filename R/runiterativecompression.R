runiterativecompression <- function(t, x, segmentation_threshold = segmentation_threshold, verbose = F){
  if(verbose){
    print("Running iterative compression to segment read depth data")
  }
  print(paste0("Initial segment count: ", nrow(t)))
  converged <- F
  compress <- t
  if(nrow(t) == 1){
    converged = T
  }
  while(!converged){
    windows <- nrow(compress)
    compress <- compressdata(compress, x, segmentation_threshold)
    if(nrow(compress) == windows | nrow(compress) == 1){
      converged <- T
    }
  }
  t$segment <- findInterval(t$pos, vec = compress$pos)
  return(t)
}

