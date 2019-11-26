callbreakpoints <- function(t, predictedpositions, maxpeak){
  t$breakpoint <- F
  for(window in 1:(nrow(t)-1)){
    if(t$segment[window] != t$segment[window+1]){
      if(t$CN[window] == t$CN[window+1]){
        next
      }
      diffs <- abs(c(min(abs(t$residual[window] - predictedpositions + maxpeak)), min(abs(t$residual[window+1] - predictedpositions + maxpeak))))
      decision <- which.max(diffs) - 1
      t$breakpoint[window + decision] <- T
    }
  }
  return(t)
}
