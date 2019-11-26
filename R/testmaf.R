#' @export
testMAF <- function(CN, tp){
  ## This special case results in division by zero
  if(CN == 0 & tp == 1){
    return(c("0" = 0.5))
  }
  #if(CN < 0 | CN > 11){
  #  stop("Please enter a CN between 0 and 8")
  #}
  np <- 1-tp
  halfcn <- ceiling(CN/2)
  if(CN < 15){
    major_allele_possibilities = seq(from = 0, to = min(10, CN), by = 1)
  }
  else if(CN < 60){
    major_allele_possibilities = c(seq(from = 0, to = 10, by = 1), seq(from = 15, to = min(50, CN), by = 5))
  }
  else if(CN < 150){
    major_allele_possibilities = c(seq(from = 0, to = 10, by = 1), seq(from = 15, to = 50, by = 5), seq(from = 60, to = min(140, CN), by = 10))
  }
  else if(CN < 250){
    major_allele_possibilities = c(seq(from = 0, to = 10, by = 1), seq(from = 15, to = 50, by = 5), seq(from = 60, to = 140, by = 10), seq(from = 150, to = CN, by = 50))
  }
  else{
    major_allele_possibilities = c(seq(from = 0, to = 10, by = 1), seq(from = 15, to = 50, by = 5), seq(from = 60, to = 140, by = 10), seq(from = 150, to = 240, by = 50), seq(from = 250, to = CN, by = 100))
  }
  output <- c()
  for(al in major_allele_possibilities){
    majoraf <- ((al * tp) + (1 * np))/((CN * tp) + (2 * np))
    output[as.character(al)] <- majoraf
  }
  return(output)
}
#' @export
testMAF_sc <- function(CN, tp){
  fraction = CN - floor(CN)
  if(fraction == 0){
    return(testMAF(CN, tp))
  }
  np = 1-tp
  base_cn <- floor(CN)
  which_fraction = c(0, 1)
  
  major_allele_possibilities = seq(from = 0, to = base_cn, by = 1)
  major_allele_possibilities = sort(c(major_allele_possibilities, major_allele_possibilities + fraction))
  output <- c()
  for(al in major_allele_possibilities){
    output[paste0(al)] <- ((al * tp) + (1 * np))/((CN * tp) + (2 * np))
  }
  return(output)
}

