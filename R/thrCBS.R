#' Detection CNA
#'
#' Detection CNA using MAD or SD.
#' 
#' @param data CBS data
#' @param method method of detection threshold. Default \code{MAD}
#' @param the numeric, if need to establishe custom threshold
#' @param n numeric, multiple to method 
#' 
#' @export


thrCBS <- function(data, method=c("MAD", "SD"), thr=NULL, n=1) {
  if(!inherits(data, "DNAcopy"))
    stop("Data must be DNAcopy")
  
  total = ncol(data$data)-2    
  res = data.frame(matrix(vector(), 0, 6, dimnames=list(c(), colnames(data$output))), 
                   stringsAsFactors=F)
  
  for(smp in 1:total){    
    forTHR = subset(data, c(1:22), smp)
    datSEG = subset(data, c(1:22,"X","Y"), smp)
    
    if(!is.null(thr)){
      zero = median(as.matrix(forTHR$data[,3]))
      x  = c(zero+thr, -thr+zero)
    } else {
      x <- thrCalc(forTHR, method=match.arg(method),n=n)
    }
    
    resTMP <- datSEG$output[which(datSEG$output$seg.mean >= x[1] | datSEG$output$seg.mean <= x[2]),]
    res = rbind(res, resTMP)
  }
  #return with exclude identical start and end
  res = res[(res$loc.end - res$loc.start) != 0,]
  xRs = list()
  for(i in unique(res$ID)){
    xRs[[i]] = res[res$ID == i,]
  }
  class(xRs) = c("list", "CBS")
  xRs
}
