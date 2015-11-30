#' Calculate threshold
#' 
#' find threshold
#' 
#' @param data subset of CBS
#' @param method perform by \code{MAD} or \code{SD}
#' @param numeric, multiple to method
#' 
#' @export



thrCalc <- function(data, method=c("MAD", "SD"), n=1){
  if(inherits(data, 'DNAcopy')){
    y <- as.matrix(na.omit(data$data[,3])) 
  } else if (inherits(data, 'numeric')) {
    y <- as.matrix(na.omit(data))
  } else {
    stop("Sorry, I can not calculate this") 
  }
  mediana <- median(y)
  method = match.arg(method)
  switch(method, 
         SD = {
           x <- c(mediana+(n*sd(y)), mediana-(n*sd(y)))
         },
         MAD = {
           x <- c(mediana+(n*mad(y)), mediana-(n*mad(y)))
         }
  )
  return(x)  
}