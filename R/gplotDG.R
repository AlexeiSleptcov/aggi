#' Plot dendrogramm for Agilent microarray data
#' 
#' Clustering tree by Hierarchical clustering algorithm.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param methodDist method for Distance Matrix Computation. See \code{\link{dist}}. Default: {"euclidean"}.
#' @param methodHclust method for hierarchical cluster analysis. See \code{\link{hclust}}. Default: {"complete"}.
#' @param rotate if TRUE, rotates plot by 90 degrees. Default: {TRUE}.
#' @param title title
#' @param ... see \code{\link{ggdendrogram}}
#'  
#' @export

gplotDendro <- function(data, 
                   methodDist = "euclidean",
                   methodHclust = "complete",
                   rotate = TRUE,
                   title = NULL,
                   ...){
  
  if(inherits(data,"RGList")){
    dannye = na.omit( suppressWarnings(log10(data$R/data$G)) )
  } else if(inherits(data,"MAList")) {
    dannye = data$M
  } else if(inherits(data,c("EList", "EListRaw"))) {
    dannye = data$E
  } else {
    stop("Data must RGList, MAList or EList class object",  call. = FALSE)
  }
  
  if(is.null(title))
    mytitle = paste("Distance:", methodDist, "/", "Hclust:", methodHclust)
  else 
    mytitle = title
  
  hc = hclust( dist( t(dannye), method = methodDist), method = methodHclust)
  ggdendrogram(hc, rotate = rotate, size = 2, hang=-1,  xlab = methodDist, ...) + 
    labs(title = mytitle)
}

