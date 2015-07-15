
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
#' @param print print result K-Means clustering
#'  
#' @export

plotKM <- function(data, 
                   centers = 2,
                   methodDist = "euclidean",
                   print = FALSE){
  
  if(inherits(data,"RGList")){
    dannye = log10(data$R/data$G)
  } else if(inherits(data,"MAList")) {
    dannye = data$M
  } else if(inherits(data,"EList")) {
    dannye = data$E
  } else {
    stop("Data must RGList, MAList or EList class object")
  }
  
  MDS.data = cmdscale( dist( t(dannye), method = methodDist));
  km = kmeans(MDS.data, centers, nstart=1000)
  plot(MDS.data, type="n")
  text(MDS.data, rownames(MDS.data), col = km$cluster)
  points(km$centers, col = 1:centers, pch = 8)
  
  if(print){
    kmr = list()
    for(i in 1:centers){
      kmr[[i]] = names(which(km$cluster == i))
    }
    return(kmr)
  }
}