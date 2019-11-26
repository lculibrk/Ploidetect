get_coverage_characteristics <- function(tp, ploidy, maxpeak){
  np=1-tp
  predicted_diff <- (maxpeak - np*maxpeak)/(ploidy - np*(ploidy - 2))
  predicted_homd <- maxpeak - predicted_diff*ploidy
  predictedpositions <- seq(from = predicted_homd, by = predicted_diff, length.out = 11)
  names(predictedpositions) <- 0:10
  return(list("maxpeak" = maxpeak,
              "ploidy" = ploidy,
              "tp" = tp,
              "homd" = predicted_homd,
              "diff" = predicted_diff,
              "cn_by_depth" = predictedpositions))
}
