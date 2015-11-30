#' Map location
#'
#' Crosome and map
#' 
#' @param data MAList or RGList
#' @param obj ligical, create \code{chrom} and \code{maploc} objects
#' @param exclude exclude probes: \code{"default"} exclude all probes besides regular probes.
#' \code{"control"} - exclude only control probes = 0 and 1.
#' \code{"probe"} - exclude your probes defined in \code{probe} argument.
#' \code{"none"} - no excluding
#' @param probe numeric, which probe you need exclude
#' @param verbose logical, system messages
#' @param include logical, if TRUE drop any probes described in non-chromosome (1:24,X,Y) regions
#' @param ingenes logical, if TRUE create \code{maploc} and \code{chrom} column in result data
#'
#' @export

chrmap <- function(data, obj = TRUE, exclude="default", probe=-15000,
         verbose=1, include=FALSE, ingenes=FALSE){
  #Check data and parameters
  if(!inherits(data, c('RGList', 'MAList')))
    stop("data must be RGList or MAList")
  if(!inherits(obj, 'logical')) #check obj
    stop("obj must be logical")
  if(!inherits(probe, 'numeric')) #check probe
    stop("probe must be numeric")
  if(!inherits(verbose, 'numeric')) #check probe
    stop("verbose must be numeric")
  
  #Create only chromosome information
  chrom <- gsub("chr([0-9XY]+):.*", "\\1", data$genes[, "SystematicName"])
  chrom <- gsub("chr1_gl000192_random:.*", "\\1", chrom) 
  #Drop probes
  if(include){
    dropProbe = 1
    drop <- which(data$genes[, "ControlType"] != probe)
  } else {
    #check probes
    up = unique(data$genes$ControlType)
    if(exclude == "control")
      if(!(1 %in% up || -1 %in% up))
        exclude = "none"
      #switch
      switch(exclude, 
             default = {
               dropProbe <- which(data$genes[, "ControlType"] != 0)
             },
             probe = {
               dropProbe = 1
               drop <- which(data$genes[, "ControlType"] == probe)
             },
             control = {
               dropProbe <- which(data$genes[, "ControlType"] %in% c(1,-1))
             },
             none = {
               dropProbe <- NULL
             },
             stop("exclude is incorrect")
      )
  } 
  
  #check dropProbe
  if(length(dropProbe)!=0) {
    #Drop probes
    if(exclude != "probe" && !include)
      drop <- c(which(!chrom %in% c(1:22, "X", "Y")), dropProbe )
    result = data[-drop, ]
    if(verbose == 1){
      cat(nrow(data$genes)-nrow(result$genes), "probes were excluded\n")
    }
    include = data$genes[-drop, "SystematicName"]
  } else {    
    #No drop
    if(exclude != 'none'){
      cat("No exclude probes\n")
    }
    include = data$genes[, "SystematicName"]
    result = data
  }
  #Create objects
  if(obj){
    chromDrop <- gsub("chr([0-9XY]+):.*", "\\1", include)
    chrom <<- ordered(chromDrop, levels=c(1:22,"X","Y"))
    maploc <<- as.numeric(gsub(".*:([0-9]+)-.*", "\\1", include))
  } else {
    if(ingenes){
      chromDrop <- gsub("chr([0-9XY]+):.*", "\\1", include)
      result$genes$chr = ordered(chromDrop, levels=c(1:22,"X","Y"))
      result$genes$maploc = as.numeric(gsub(".*:([0-9]+)-.*", "\\1", include))
    }
  }
  result
}
