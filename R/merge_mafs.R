merge_mafs <- function(inputmafs, na.rm = T, flip = F, exp = F){
  if(exp){
    if(length(na.omit(inputmafs)) == 0){
      return(NA_character_)
    }
    return(paste(na.omit(inputmafs), collapse = ";"))
  }
  if(na.rm){
    inputmafs <- inputmafs[which(!is.na(inputmafs))]
  }
  ## Compute magnitude about 0.5
  mafs_about_zero <- inputmafs - 0.5
  negatives <- which(mafs_about_zero < 0)
  positives <- which(mafs_about_zero > 0)
  ## Compute absolute magnitude about 0.5 and add 0.5
  flipped_mafs <- abs(mafs_about_zero)
  if(length(positives) > length(negatives)){
    out <- mean(0.5 + flipped_mafs)
  }else{
    out <- mean(0.5 - flipped_mafs)
  }
  if(flip){
    out <- abs(out - 0.5) + 0.5
  }
  return(out)
}
unmerge_mafs <- function(merged_mafs, flip = F){
  mafs <- na.omit(merged_mafs)
  if(flip){
    return(abs(as.numeric(unlist(lapply(mafs, strsplit, split = ";"), recursive = T)) - 0.5) + 0.5)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}
flip_merged_mafs <- function(merged_mafs){
  mafs <- unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
  flipped_mafs <- lapply(mafs, function(x)merge_mafs(abs(as.numeric(x) - 0.5) + 0.5, exp = T))
  return(flipped_mafs)
}

unmerge_mafs_grouped <- function(merged_mafs, flip = F){
  if(flip){
    mafs <- unlist(lapply(merged_mafs, strsplit, split = ";"), recursive = F)
    mafs <- lapply(mafs, function(x){median(abs(as.numeric(x) - 0.5) + 0.5)})
    return(mafs)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}
