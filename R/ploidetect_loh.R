ploidetect_loh <- function(purity, CNA_object){
  CNA_object$CN <- round(CNA_object$CN)
  CNA_object_mafs <- CNA_object[which(!is.na(CNA_object$maf)),]
  CNA_object_nomafs <- CNA_object[which(is.na(CNA_object$maf)),]
  if(nrow(CNA_object_nomafs) > 1){
    CNA_object_nomafs$LOH <- F
  }

  allCNs <- unique(CNA_object_mafs$CN)
  allCNs <- allCNs[allCNs > 1 & allCNs <= 8]
  maf_possibilities <- list()
  for(copynumber in allCNs){
    mafs <- testMAF(copynumber, tp = purity)
    LOH_maf <- mafs[length(mafs)]
    almost_LOH_maf <- mafs[length(mafs) - 1]
    maf_possibilities[paste0(copynumber)] <- list(c(LOH_maf, almost_LOH_maf))
  }
  segments <- CNA_object_mafs %>% group_by(chr, segment) %>% dplyr::summarise(mafs = merge_mafs(maf, exp = T), median_segment = median(median_segment), CN = median(CN))
  
  segments$LOH <- F
  for(segment in 1:nrow(segments)){
    CN <- paste0(segments$CN[segment])
    if(CN == "1"){
      segments$LOH[segment] <- T
    }
    if(CN %in% names(maf_possibilities)){
      candidates <- unlist(maf_possibilities[CN])
      decision <- which.min(abs(median(unmerge_mafs(segments$mafs[segment], flip = T)) - candidates))
      if(length(decision) == 0){
        segments$LOH[segment] <- F
      }else if(decision == 1){
        segments$LOH[segment] <- T
      }else{
        segments$LOH[segment] <- F
      }
    }
  }
  
  CNA_object_mafs <- left_join(CNA_object_mafs, segments[,c("chr", "segment", "LOH")], by = c("chr", "segment"))
  if(nrow(CNA_object_nomafs) > 0){
    CNA_object <- rbind.data.frame(CNA_object_mafs, CNA_object_nomafs) %>% arrange(chr, pos)
  }else{
    CNA_object <- CNA_object_mafs %>% arrange(chr, pos)
  }
  state_call <- CNA_object %>% ungroup %>% group_by(chr, segment) %>% dplyr::summarise(CN = first(CN), LOH = first(LOH))
  ## Assign HMM-like states to CN/LOH combinations
  state_call$state <- 0
  for(segment in 1:nrow(state_call)){
    ## 0 = HOMD, 1 = 1cp, 2 = 2cp het, 3 = 2cp HOM, 4 = 3cp het, 5 = 3cp HOM, 6 = 4cp het, 7 = 4cp HOM, 8 = 5cp+
    CN <- state_call$CN[segment]
    if(CN >= 5){
      state_call$state[segment] <- 8
    }else if(CN == 0){
      state_call$state[segment] <- 0
    }else if(CN == 1){
      state_call$state[segment] <- 1
    }else if(CN == 2){
      if(state_call$LOH[segment]){
        state_call$state[segment] <- 3
      }else{
        state_call$state[segment] <- 2
      }
    }else if(CN == 3){
      if(state_call$LOH[segment]){
        state_call$state[segment] <- 5
      }else{
        state_call$state[segment] <- 4
      }
    }else if(CN == 4){
      if(state_call$LOH[segment]){
        state_call$state[segment] <- 7
      }else{
        state_call$state[segment] <- 6
      }
    }
  }
  CNA_object <- CNA_object %>% ungroup %>% left_join(state_call[,c("chr", "segment", "state")], by = c("chr", "segment"))
  return(CNA_object)
}
