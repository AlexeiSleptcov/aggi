#' Plot dendrogramm for Agilent microarray data
#' 
#' Clustering tree by Hierarchical clustering algorithm.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param cyanine only for RGList data class.
#' 
#' @param methodDist method for Distance Matrix Computation. See \code{\link{dist}}. Default: {"euclidean"}.
#' @param methodHclust method for hierarchical cluster analysis. See \code{\link{hclust}}. Default: {"complete"}.
#'  
#' @export

plotDG <- function(data, 
                   cyanine = c("R","G"),
                   methodDist = "euclidean",
                   methodHclust = "complete"){
  
  if(inherits(data,"RGlist")){
    cyanine = match.arg(cyanine)
    dannye = data[[cyanine]]
  } else if(inherits(data,"MAlist")) {
    dannye = data$M
  } else if(inherits(data,"Elist")) {
    dannye = data$E
  } else {
    stop("Data must RGList, MAList or EList class object")
  }
  
  if(require("limma"))
    plot(hclust( dist( t(dannye), method = methodDist), method = methodHclust), hang=-1)
}


#' Plot clustering MDS for Agilent microarray data
#' 
#' Clustering MDS by K-Means clustering algorithm.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param cyanine only for RGList data class.
#' @param centers either the number of clusters, say \emph{k}, or a set initial (distinct)
#' cluster centers. If a number, a random set of (distinct) rows in x is chosen as initial centres.
#' See \code{\link{kmeans}}. Default: \code{2}.
#' 
#' @param methodDist method for Distance Matrix Computation. See \code{\link{dist}}. Default: {"euclidean"}.
#' @param methodHclust method for hierarchical cluster analysis. See \code{\link{hclust}}. Default: {"complete"}.
#' 
#' @param print print result K-Means clustering
#'  
#' @export

plotDG <- function(data, 
                   cyanine = c("R","G"),
                   centers = 2,
                   methodDist = "euclidean",
                   methodHclust = "complete",
                   print = FALSE){
  
  if(inherits(data,"RGlist")){
    cyanine = match.arg(cyanine)
    dannye = data[[cyanine]]
  } else if(inherits(data,"MAlist")) {
    dannye = data$M
  } else if(inherits(data,"Elist")) {
    dannye = data$E
  } else {
    stop("Data must RGList, MAList or EList class object")
  }
  
  MDS.data = cmdscale( dist( t(dannye), method = methodDist), method = methodHclust);
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