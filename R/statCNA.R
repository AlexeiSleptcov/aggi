#' simple CNA statistic
#' 
#' @param data CBS data
#' 
#' @export 

statCNA <- function(data){
  if(inherits(data,"DNAcopy"))
    data = thrCBS(data)
  if(!inherits(data,"CBS"))
    stop("First, data must be CBS or DNAcopy", call. = 0)
  
  res = list()
  for(i in names(data)){
    data[[i]]$chrom = ordered(data[[i]]$chrom, c(1:22,"X","Y"))
    data[[i]] = list(data[[i]][which(data[[i]]$seg.mean > 0),],
                     data[[i]][which(data[[i]]$seg.mean < 0),])
    names(data[[i]]) = c("gain", "loss")
    res[[i]] = t(data.frame(gain = c(summary(data[[i]]$gain$chrom), length(data[[i]]$gain$chrom)), 
                            loss = c(summary(data[[i]]$loss$chrom), length(data[[i]]$loss$chrom))))
  }
  res
}