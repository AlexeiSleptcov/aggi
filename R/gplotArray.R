#' Plot array
#' 
#' @param data RGList, MAList, EListRaw or EList
#' @param smp array
#' 
#' @export


gplotArray <- function(data, smp = 1){
  if(inherits(data,"RGList")){
    x = MA.RG(data)$M[,smp]
  } else if (inherits(data, "MAList")){
    x = data$M[,smp]
  } else if (inherits(data,c("EListRaw","EList"))){
    x = data$E[,smp]
  } else stop("First, data must be RGList or EListRaw", call.=FALSE)
  
  b = data.frame(row = data$genes$Row, col = data$genes$Col, log = x)
  
  if( inherits(data,c("EListRaw","EList")) ){
    colours = c("#DAF07B","#6ABE54","#066730")
  } else colours = c("#056730","#FFFDB5","#A61524")
  
  ggplot(b, aes(x=row, y=col, fill=log)) + geom_raster() + 
    ggtitle( colnames(data)[smp] ) + 
    scale_fill_gradientn(colours=colours, na.value = "black") +
    theme(
      panel.grid.major = element_line(colour = "#252525"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "black")
      ) + 
    theme(legend.position="none")
}