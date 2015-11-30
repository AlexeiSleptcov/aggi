#' Derivative log ratio spread
#'
#' Derivative log ratio spread (DLRS) 
#' is a measurement of point-to-point consistency or noisiness in log ratio data
#' 
#' @param data RGList or MAList
#' 
#' @return Excellent < 0.2, good < 0.3, bad > 0.3
#' 
#' @export

DLRSpread <- function(data) {
  #check data 
  if(inherits(data, "RGList")){
    presto = MA.RG(data)
  } else if (inherits(data, "MAList")){
    presto = data
  } else {
    stop("Data must be MAList or RGList class")
  }
  #create objects
  proto = chrmap(presto, obj=T,exclude='control')
  logratio = proto$M
  arrayNames = colnames(logratio)
  #create object for DLRS
  vecs = cbind(chrom,maploc,logratio)
  colnames(vecs) = c("chr","pos",arrayNames)
  vecpos = vecs[order(vecs[,2]),]
  vcp = vecpos[order(vecpos[,1]),]  
  # Run this on the Log ratios (ordered by chromosome and position).
  dlrs = matrix(NA,nrow=2,ncol=length(arrayNames),dimnames=list(c('DLRSpread',"Quality"),arrayNames))
  for(i in 3:ncol(vcp)){
    nx <- length(vcp[,i])
    if (nx<3) {
      stop("Vector length>2 needed for computation")
    }
    tmp <- embed(vcp[,i],2)
    diffs <- as.numeric(na.omit(tmp[,2]-tmp[,1]))
    assess = IQR(diffs)/(sqrt(2)*1.34)
    dlrs[1,i-2] <- round(assess, 3)
    dlrs[2,i-2] <- if(assess < 0.2 ) {
      "Excelent"
    } else if (assess < 0.3) {
      "Good"
    } else {
      'Bad'
    }
  }
  dlrs
}