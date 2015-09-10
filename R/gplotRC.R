#' Plot Row&Column
#'
#' Plot of median distribution by rows and columns
#' 
#' @param data RGarray, Earray
#' @param smp array
#' 
#' @export

gplotRC  <- function(data, smp=1){
  
  #http://stackoverflow.com/questions/17370460/scatterplot-with-alpha-transparent-histograms-in-r
  theme0 <- function(...) theme( legend.position = "none",
                                 panel.margin = unit(0,"null"),
                                 axis.ticks = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.ticks.length = unit(0,"null"),
                                 axis.ticks.margin = unit(0,"null"),...)
  if(!inherits(data, c("RGList", "EListRaw", "EList")))
    stop("Data must be RGList or EListRaw")
  data = makeArrays(data)
  
  if(inherits(data, "RGarray")){
    
    negR = R = data$Arrays[,,smp, 1]
    negG = G = data$Arrays[,,smp, 2]
    
    negR[data$ControlType != -1] = NA
    negG[data$ControlType != -1] = NA
    
    dCol = data.frame(
      Col         = c(1:dim(R)[2], 
                      1:dim(G)[2],
                      1:dim(R)[2],
                      1:dim(G)[2]
                      ),
      Intensities = c(apply(R, 2, median, na.rm = TRUE), 
                      apply(G, 2, median, na.rm = TRUE),
                      apply(negR, 2, median, na.rm = TRUE), 
                      apply(negG, 2, median, na.rm = TRUE)
                      ),
      Probes        = c(rep("Regular R", dim(R)[2]), 
                        rep("Regular G", dim(G)[2]),
                        rep("Negatives R", dim(R)[2]), 
                        rep("Negatives G", dim(G)[2])
                      )
    )
    
    levels(dCol$Probes) = c("Regular R","Regular G" ,"Negatives R","Negatives G")
    
    dRow = data.frame(
      Row         = c(1:dim(R)[1], 
                      1:dim(G)[1],
                      1:dim(R)[1],
                      1:dim(G)[1]
      ),
      Intensities = c(apply(R, 1, median, na.rm = TRUE), 
                      apply(G, 1, median, na.rm = TRUE),
                      apply(negR, 1, median, na.rm = TRUE), 
                      apply(negG, 1, median, na.rm = TRUE)
      ),
      Probes        = c(rep("Regular R", dim(R)[1]), 
                        rep("Regular G", dim(G)[1]),
                        rep("Negatives R", dim(R)[1]), 
                        rep("Negatives G", dim(G)[1])
      )
    )
    
    levels(dRow$Probes) = c("Regular R","Regular G" ,"Negatives R","Negatives G")
    
    top <- ggplot(dCol, aes(x=Col, y=Intensities, fill=Probes)) + 
      stat_smooth(method = loess) + 
      scale_fill_manual(values = c("#066730","#A61625", "#066730","#A61625")) +
      scale_colour_manual(values = c("#066730","#A61625", "#066730","#A61625")) +
      geom_point(aes(colour = Probes)) +
      ggtitle(paste("CGH Array N", smp, sep="")) 
    
    bottom <- ggplot(dRow, aes(x=Row, y=Intensities, fill=Probes)) + 
      stat_smooth(method = loess) + 
      scale_fill_manual(values = c("#066730","#A61625", "#066730","#A61625")) +
      scale_colour_manual(values = c("#066730","#A61625", "#066730","#A61625")) +
      geom_point(aes(colour = Probes))
    
    suppressWarnings(grid.arrange(arrangeGrob(top), arrangeGrob(bottom)))
    
  } else if (inherits(data, "Earray")){
    
    neg = x = data$Arrays[,,smp]
    neg[data$ControlType != -1] = NA
    
    dCol = data.frame(
      Col         = c(1:dim(x)[2], 
                      1:dim(neg)[2]),
      Intensities = c(apply(x, 2, median, na.rm = TRUE),
                      apply(neg, 2, median, na.rm = TRUE)),
      Probes        = c(rep("Regular", dim(x)[2]),
                        rep("Negatives", dim(neg)[2]))
    )
    
    levels(dCol$Probes) = c("Regular", "Negatives")
    
    dRow = data.frame(
      Row         = c(1:dim(x)[1], 
                      1:dim(neg)[1]),
      Intensities = c(apply(x, 1, median, na.rm = TRUE),
                      apply(neg, 1, median, na.rm = TRUE)),
      Probes        = c(rep("Regular", dim(x)[1]),
                        rep("Negatives", dim(neg)[1]))
    )
    
    levels(dRow$Probes) = c("Regular", "Negatives")
    
    top <- ggplot(dCol, aes(Col, Intensities, colour=Probes)) + 
      stat_smooth(method = loess, fill="#066730") + 
      geom_point() +
      ggtitle(paste("Expression Array N", smp, sep="")) 
    
    bottom <- ggplot(dRow, aes(Row, Intensities, colour=Probes)) + 
      stat_smooth(method = loess, fill="#066730") + 
      geom_point()
    
    suppressWarnings(grid.arrange(arrangeGrob(top), arrangeGrob(bottom)))
    
  } else stop("Data musr RGarray or Earray", call. = FALSE)
}
