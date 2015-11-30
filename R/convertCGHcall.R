#' Convert MAlist to CGHcall
#'
#' Convert MAlist to CGHcall
#'
#' @param x MAList
#' @param XY logical, if TRUE exclude XY
#' @param numChr number to chr
#'
#' @export

convertCGHcall <- function(x, XY = T, numChr= T){
  
  if(!inherits(x, "MAList"))
    stop("First, data must be  MAList", call. = FALSE)
  
  x = chrmap(x, obj = FALSE, verbose = 0, ingenes = FALSE)
  x = x[-anyDuplicated(x$genes$SystematicName),]
  if(numChr){
    x$genes$Chr<- lapply(strsplit(as.character(x$genes$SystematicName), "\\chr"), "[", 2)
    x$genes$Chr<- lapply(strsplit(as.character(x$genes$Chr), "\\:"), "[", 1) 
  } else {
    x$genes$Chr<- lapply(strsplit(as.character(x$genes$SystematicName), "\\:"), "[", 1)
  }
  x$genes$Position <- lapply(strsplit(as.character(x$genes$SystematicName), "\\:"), "[", 2)
  x$genes$Position <- lapply(strsplit(as.character(x$genes$Position), "\\-"), "[", 1)
  x$genes$Chr <- unlist(x$genes$Chr)
  x$genes$Position <- unlist(x$genes$Position)
  x$genes$Chr[x$genes$Chr == '1_gl000192_random'] <- NA
  if(XY){
    x$genes$Chr[x$genes$Chr == 'X'] <- 23
    x$genes$Chr[x$genes$Chr == 'Y'] <- 24
    x$genes$Chr <- as.integer(x$genes$Chr)    
  }
  x$genes$Position <- as.integer(x$genes$Position)
  y <- data.frame("PROBE" = x$genes$ProbeName,
             "CHROMOSOME" = x$genes$Chr,
             "START_POS" = x$genes$Position,
             "STOP_POS" = x$genes$Position+60) 
  cbind(y,x$M)
}