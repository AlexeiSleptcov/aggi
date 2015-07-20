#' MA-plot: plot differences versus averages for Agilent microarray data
#' 
#' A generic function which produces an MA-plot for an object containing microarray data.
#' 
#' @param data microarray data. Only for \code{RGlist}, \code{MAlist} or \code{ELits} class objects.
#' @param array number array or name array to be displayed. Default \code{1}.
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
#' 
plotMA.agi <- function(data, 
                       array = 1, 
                       probe = NULL,
                       alpha = 0.5,
                       size = 2){
  if(!inherits(data, c("RGList", "MAList")))
    stop("First, data must be RGList or MAList", call. = FALSE)
  
  if(inherits(data, "RGList"))
    data = limma::MA.RG(data)
  
  # predata
  y.lim  = c(min(data$M[,array], na.rm = 1), max(data$M[,array], na.rm = 1))
  x.lim  = c(min(data$A[,array], na.rm = 1), max(data$A[,array], na.rm = 1))
  if(is.null(probe))
    data = na.omit(cbind(data$M[,array], data$A[,array], data$genes$ControlType))
  else {
    whi  = which(data$genes$ControlType == probe)
    data = na.omit(cbind(data$M[whi,array], data$A[whi,array], data$genes$ControlType[whi]))
  }
  
  colnames(data) = c("M", "A", "ProbeType")
  
  # ggplot
  ggplot2::ggplot(as.data.frame(data), aes(A,M)) +
    geom_point(aes(color = factor(ProbeType)), alpha = alpha, size = size) +
    guides(color = guide_legend(title = "Probe type", override.aes = list(size = 4)),
           alpha = guide_legend(label = 0)) +
    labs(title = colnames(data$M)[array])
}




