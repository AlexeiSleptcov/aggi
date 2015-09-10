#' Remove Outliners
#' 
#' Detect and remove outliners by intra-array or inter-array methods.
#'
#' @param data RGList, EListRaw
#' @param flags logical, if \code{TRUE} find outliners by flags else by ratio of mean and median
#' @param method if \code{intra} - detect and remove outliners intra-array, if \code{between} - 
#' detect and remove outliners inter-array.
#' 
#' @export


rmOutliners <- function(data, flags = TRUE,
                        method=c("intra", "between")){
  if(!inherits(data,c("RGList", "EListRaw")))
    stop("First, data must be RGList or EListRaw", call.=FALSE)
  
  method = match.arg(method)
  if(inherits(data, "RGList")){
    a = c()
    
    if(method == "intra"){
      for(i in 1:ncol(data)){
        a = findOutliners(data, smp = i, flags)
        data$R[a,i] = NA
        data$G[a,i] = NA
      }
    } else {
      for(i in 1:ncol(data)){
        a = c(a, findOutliners(data, smp = i, flags))
      }
      data$R[unique(a),] = NA
      data$G[unique(a),] = NA
    }
  } else {
    a = c()
    
    if(method == "intra"){
      for(i in 1:ncol(data)){
        a = findOutliners(data, smp = i, flags)
        data$E[a,i] = NA
      }
    } else {
      for(i in 1:ncol(data)){
        a = c(a, findOutliners(data, smp = i, flags))
      }
      data$E[unique(a),] = NA
    }
  }
  data
}