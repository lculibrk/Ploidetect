#' Compute coverage characteristics from tp/ploidy/deptht
#' 
#' \code{get_coverage_characteristics} uses tp/ploidy/depth information to model expected shifts in depth
#' per copy number, predicted HOMD depth, and returns these estimates with predictions of depth for integer
#' copy numbers from 0 to 10. Used in Ploidetect as a function to mainly return this vector.
#' 
#' @param tp tumour purity (float)
#' @param ploidy ploidy estimate (integer)
#' @param maxpeak most common read depth in the tumour - found using kernal density estimation
#' @return a named list of variables
get_coverage_characteristics <- function(tp, ploidy, maxpeak){
  ## Normal tumor purity
  np=1-tp
  ## Formula for depth differential from 1-copy change
  ## I have discovered a truly marvellous derivation for this, which
  ## this comment section is too small to contain
  predicted_diff <- (maxpeak - np*maxpeak)/(ploidy - np*(ploidy - 2))
  ## Predict the homozygous deletion depth
  predicted_homd <- maxpeak - predicted_diff*ploidy
  ## Create a vector of read depth by copy number
  predictedpositions <- seq(from = predicted_homd, by = predicted_diff, length.out = 11)
  names(predictedpositions) <- 0:10
  ## Return list
  return(list("maxpeak" = maxpeak,
              "ploidy" = ploidy,
              "tp" = tp,
              "homd" = predicted_homd,
              "diff" = predicted_diff,
              "cn_by_depth" = predictedpositions))
}
