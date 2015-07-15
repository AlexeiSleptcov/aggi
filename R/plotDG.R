#' Plot dendrogramm for Agilent microarray data
#' 
#' Clustering tree by Hierarchical clustering algorithm.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param methodDist method for Distance Matrix Computation. See \code{\link{dist}}. Default: {"euclidean"}.
#' @param methodHclust method for hierarchical cluster analysis. See \code{\link{hclust}}. Default: {"complete"}.
#'  
#' @export

plotDG <- function(data, 
                   methodDist = "euclidean",
                   methodHclust = "complete"){
  
  if(inherits(data,"RGlist")){
    dannye = log10(data$R/data$G)
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

