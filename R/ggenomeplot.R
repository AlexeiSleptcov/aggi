#' Genome plot 
#' 
#' Genome plot
#' 
#' @param data MAList or DNAcopy class data
#' @param smp numeric, number of sample. if NULL all sample will be saved in current directory
#' @param chr 1:23, X, Y. If NULL karyotype will be made
#' @param span numeric, the parameter α which controls the degree of smoothing.
#' @param ... params for \code{\link{qplot}} only if param \code{chr} was applied 
#' 
#' @export


ggenomeplot <- function(data, 
                       smp= NULL, 
                       chr=NULL, 
                       span=0.1,
                       ...){
  
  if(inherits(data, "MAList")){
    data = data[order(data$genes$maploc),]
    data = data[order(data$genes$chr),]
    nameArray = colnames(data)
    locData = data$M
    chrData = data$genes$chr
    mapData = data$genes$maploc
  } else if(inherits(data, "DNAcopy")) {
    nameArray = colnames(data$data)[-(1:2)]
    locData = as.data.frame(data$data[,-(1:2)])
    chrData = factor(data$data$chrom, levels=c(1:23,"X","Y"))
    mapData = data$data$maploc
  } else {
    stop("Must be DNAcopy or MAList")
  }
  
  if(!is.null(smp))
    limarg(smp, maxi=data.frame(array = length(nameArray)), maxl=1)
  
  if(!is.null(chr)){
    if(!(chr %in% c(1:22,'X','Y')))
      stop("argument chr must be from 1 to 22, X or Y")
    locData = locData[ chrData == chr,]
    xlabName = paste(chr,"Chromosome")
  } else {
    xlabName = "Genome"
  }
  
  if(!is.null(smp)){
    nam = as.double(smooth(na.omit(locData[,smp])))
    if(!is.null(chr)){
      qplot(1:length(nam), nam, geom=c("point"), ylab = nameArray[smp], xlab = xlabName, alpha = I(.2), size = I(1), ylim=c(-2,2),na.rm=TRUE, ...) +
        stat_smooth(method="loess", size = 1,span = span, level =0.99)
    } else {
      require(RColorBrewer)
      mypal = brewer.pal(9,"Set1")[1:2]
      nam = data.frame(x = as.double(smooth(na.omit(locData[,smp]))), 
                       cyl = chrData[!is.na(locData[,smp])], 
                       pos = mapData[!is.na(locData[,smp])])
      ggplot(nam, aes(pos, x, colour=factor(x>0))) + geom_point(na.rm=TRUE, size = 1) + 
        facet_wrap(~cyl, scales = "free", shrink=F) +  ylim(-2,2) + expand_limits(x=0) +
        scale_colour_manual(values=mypal) + theme(legend.position="none") +  labs(x = nameArray[smp], y = "") +
        geom_smooth(aes(pos, x), data=nam, method="loess", size = 1,span = span, level =0.99, na.rm=TRUE, colour = "black")       
      
    }
  } else {
    require(RColorBrewer)
    mypal = brewer.pal(9,"Set1")[1:2]
    for(i in 1:length(nameArray)){
      nam = data.frame(x = as.double(smooth(na.omit(locData[,i]))), 
                       cyl = chrData[!is.na(locData[,i])], 
                       pos = mapData[!is.na(locData[,i])])
      ggplot(nam, aes(pos, x, colour=factor(x>0))) + geom_point(na.rm=TRUE, size = 1) + 
        facet_wrap(~cyl, scales = "free", shrink=F) +  ylim(-2,2) + expand_limits(x=0) +
        scale_colour_manual(values=mypal) + theme(legend.position="none") +  labs(x = nameArray[i], y = "") +
        geom_smooth(aes(pos, x), data=nam, method="loess", size = 1,span = span, level =0.99, na.rm=TRUE, colour = "black")       
      ggsave(file=paste(xlabName,"_",nameArray[i],".tiff", sep=""), scale=2, type="cairo", compression="lzw")
    }
  }
}