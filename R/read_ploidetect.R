#' @export
read_ploidetect <- function(file){
  out <- read.table(file, header = F, stringsAsFactors = F)
  return(out)
}
