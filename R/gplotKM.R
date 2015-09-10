#' Plot clustering MDS for Agilent microarray data
#' 
#' Clustering MDS by K-Means clustering algorithm.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param centers either the number of clusters, say \emph{k}, or a set initial (distinct)
#' cluster centers. If a number, a random set of (distinct) rows in x is chosen as initial centres.
#' See \code{\link{kmeans}}. Default: \code{2}.
#' 
#' @param methodDist method for Distance Matrix Computation. See \code{\link{dist}}. Default: {"euclidean"}.
#' @param groups forms by groups. Only if argument \code{old} is \code{FALSE}.
#' @param getresult print result K-Means clustering
#' 
#' @param old if TRUE draw by plot, else ggplot. Default:{FALSE}
#'  
#' @export

gplotKMeans <- function(data, 
                   centers = 2,
                   methodDist = "euclidean",
                   groups = NULL,
                   getresult = FALSE,
                   old = FALSE,
                   ...){
  
  if(inherits(data,"RGList")){
    dannye = na.omit( suppressWarnings(log10(data$R/data$G)) )
  } else if(inherits(data,"MAList")) {
    dannye = data$M
  } else if(inherits(data,c("EList", "EListRaw"))) {
    dannye = data$E
  } else {
    stop("Data must RGList, MAList or EList class object")
  }
  
  MDS.data = cmdscale( dist( t(dannye), method = methodDist));
  km = kmeans(MDS.data, centers, nstart=1000)
  
  if(getresult){
    kmr = list()
    for(i in 1:centers){
      kmr[[i]] = names(which(km$cluster == i))
    }
    return(kmr)
  }
  
  if(old){
    plot(MDS.data, type="n", xlab = "Dimension 1", ylab = "Dimension 2", ...)
    text(MDS.data, rownames(MDS.data), col = km$cluster)
    points(km$centers, col = 1:centers, pch = 8)
  } else {
    if(is.null(groups)){
      plot.dat <- as.data.frame(cbind(MDS.data, cluster = km$cluster))
      plot.dat$cluster = factor(plot.dat$cluster)
      ggplot(plot.dat,  aes(V1, V2, color = cluster, label=rownames(plot.dat))) +
        geom_text() + labs(title = "MDS plot")
    } else {
      if(!is.vector(groups))
        stop("Argument groups must be vector", call. = FALSE)
      if(length(groups) != nrow(MDS.data))
        stop("Argument groups must be equal to number of arrays", call. = FALSE)
      
      plot.dat <- as.data.frame(cbind(MDS.data, cluster = km$cluster))
      plot.dat$cluster = factor(plot.dat$cluster)
      
      if(centers != 1){
        ggplot(plot.dat,  aes(V1, V2, color = cluster, shape=groups)) +
          geom_point(size=4) + labs(title = "MDS plot")
      } else {
        ggplot(plot.dat,  aes(V1, V2, shape=groups)) +
          geom_point(size=4) + labs(title = "MDS plot")
      }
    }
  }
}