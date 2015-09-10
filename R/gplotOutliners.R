#' Plot outliners in array
#' 
#' @param data RGList, EListRaw
#' @param smp array
#' @param flags logical, if \code{TRUE} find outliners by flags else by ratio of mean and median
#'
#' @export


gplotOutliners <- function(data, smp = 1, flags = TRUE){
  if(!inherits(data,c("RGList","EListRaw")))
    stop("First, data must be RGList or EListRaw", call.=FALSE)
  
  
  a = data[findOutliners(data, smp, flags),]
  x = 1:nrow(a)
  
  if(flags){
    if(inherits(data,"RGList")){
      x[a$other$rIsFeatNonUnifOL[,smp] == TRUE | a$other$gIsFeatNonUnifOL[,smp] == TRUE] = "IsFeatNonUnifOL"
      x[a$other$rIsFeatPopnOL[,smp] == TRUE | a$other$gIsFeatPopnOL[,smp] == TRUE] = "IsFeatPopnOL"
      x[a$other$rIsSaturated[,smp] == TRUE | a$other$gIsSaturated[,smp] == TRUE] = "IsSaturated"
    } else {
      x[ a$other$gIsFeatNonUnifOL[,smp] == TRUE] = "IsFeatNonUnifOL"
      x[ a$other$gIsFeatPopnOL[,smp] == TRUE] = "IsFeatPopnOL"
      x[ a$other$gIsSaturated[,smp] == TRUE] = "IsSaturated"
    }
  } else {
    x = "Bad signals"
  }
  
  
  
  b = data.frame(row = a$genes$Row, col = a$genes$Col, log = x)
  
  ggplot(b, aes(x=row, y=col, colour=factor(log))) + geom_point() + 
    ggtitle( colnames(data)[smp] ) + theme(legend.position = "bottom") +
    guides(colour=guide_legend(title="Legend"))
}