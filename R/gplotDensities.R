#' Plot Arrays Densities
#' Plot the density of values for multiple arrays on the same plot.
#'
#' @param data ElistRaw, EList, MAList or RGList
#' @param type Type of plot: \code{density}, \code{violin}, \code{box}
#' 
#' @export

gplotDensities <- function(data, 
                           type=c("density", "violin", "box")){
  
  if( !inherits(data, c("MAList","RGList", "EList", "EListRaw")) )
    stop("First, data must be RGList, MAList or EList", call.=FALSE)
  
  if( inherits(data, "RGList") ){
    if(is.null(data$R) || length(data$G) == 0)
      stop("RG is null or length zero", call.=FALSE)
    dataR = melt(log10(data$R))
    dataR$colour = paste(dataR$Var2, "R", sep="_")
    dataG = melt(log10(data$G))
    dataG$colour = paste(dataG$Var2, "G", sep="_")
    data=na.omit(rbind(dataR,dataG))
  } else if ( inherits(data,"MAList") ){
    if(is.null(data$M) || length(data$M) == 0)
      stop("Mvalue is null or length zero", call.=FALSE)
    dataM = melt(data$M)
    dataM$colour = paste(dataM$Var2, "LR", sep="_")
    data=na.omit(dataM)
  } else {
    if(is.null(data$E) || length(data$E) == 0)
      stop("Evalue is null or length zero", call.=FALSE)
    dataE = melt(log10(data$E))
    dataE$colour = paste(dataE$Var2, "E", sep="_")
    data=na.omit(dataE)
  }
  
  type=match.arg(type)
  switch(type, 
         "box" = {
           p<-ggplot(data,aes(x=colour, y=value, fill=colour)) + 
             geom_boxplot()
         },
         "density" = {
           p<-ggplot(data, aes(x=value, colour=colour)) + 
             geom_density(size=1) + 
             labs(title = "RawData")
         },
         "volcano" = {
           p<-ggplot(data, aes(x = value)) +
             stat_density(aes(ymax = ..density..,  ymin = -..density.., 
                              fill = colour, colour = colour),
                          geom = "ribbon", position = "identity") +
             facet_grid(. ~ colour) + 
             coord_flip()
         })
  
  print(p)
  
}