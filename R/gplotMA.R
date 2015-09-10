#' MA-plot: plot differences versus averages for Agilent microarray data
#' 
#' A generic function which produces an MA-plot for an object containing microarray data.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param array number array or name array to be displayed. Default \code{1}.
#' @param columnProbe column where contain probe number
#' @param probe agilent type of probe. See in "ControlType" column. Default \code{NULL}, i.e. all.
#' @param alpha transparent points (min = 1, max = 0). Default \code{0.25}.
#' @param size size of points. Default \code{2}.
#' 
#' @details An MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages 
#' (A-values). If MA is an RGList or MAList then this function produces an ordinary within-array 
#' MA-plot. If MA is an MArrayLM object, then the plot is an fitted model MA-plot in which the 
#' estimated coefficient is on the y-axis and the average A-value is on the x-axis. 
#' 
#' @return A plot is created on the current graphics device.
#' 
#' @export
gplotMA <- function(data, 
                     array = 1, 
                     columnProbe = "ControlType",
                     probe = NULL,
                     alpha = 0.5,
                     size = 2){
  if(!inherits(data, c("RGList", "MAList")))
    stop("First, data must be RGList, MAList", call. = FALSE)
  
  if(inherits(data, "RGList"))
    data = MA.RG(data)
  
  if( is.null(ncol( data$M )) ){
    M = data[["M"]]
    A = data[["A"]]
    dataNames = "MA plot"
  } else {
    M = data[["M"]][,array]
    A = data[["A"]][,array]
    dataNames = colnames(data[["M"]])[array]
  }
  
  if(is.null(data$genes[[columnProbe]])){
    attr(MA$genes$Status, "values")
    attr(MA$genes$Status, "col")
    
    
    columnProbe = "LogIntM"
    probes = rep("medium", length(M))
    probes[M < summary(M)[2]] = "low"
    probes[M > summary(M)[4]] = "high"
    probes = factor(probes, levels = c("high", "medium", "low"))
  } else {
    if(is.null(probe))
      probes = data$genes[[columnProbe]]
    else {
      probes = data$genes[ which(data$genes[[columnProbe]] == probe) , columnProbe ]
    }
  }
  
  data = na.omit(data.frame(M, A, probes, stringsAsFactors = FALSE))
  # ggplot
  ggplot(data,  aes(A, M)) +
    geom_point(aes(colour = factor(probes) ), alpha = alpha, size = size) +
    guides(color = guide_legend(title = columnProbe, override.aes = list(size = 4)),
           alpha = guide_legend(label = 0)) +
    labs(title = dataNames) + annotation_logticks(alpha = 0.5) 
}




