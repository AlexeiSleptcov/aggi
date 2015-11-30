#' Limarg
#' 
#' Internal function. Limits of arguments
#' 
#' @param arg numeric or vector, argument
#' @param maxi check argument with maximums
#' @param maxl maxwell
#' @param class check classes
#' 
#' @export
limarg <- function(arg=NULL, maxi=NULL, maxl=NULL, class=c("numeric","integer","double")){
  if(is.null(arg))
    stop(paste("argument", substitute(arg) ,"is null"))
  if(!any(class %in% class(arg)))
    stop(paste("\nargument", substitute(arg) ,"is not", class))
  if(!is.null(maxi))
    for(i in 1:length(maxi)){
      if(arg > maxi[i]){
        nameMax = maxi[i]
        if(!is.null(names(maxi[i])))
          nameMax = names(maxi[i])
        stop(paste("argument", substitute(arg) ,"is greater then", nameMax))
      }
    }
  if(!is.null(maxl))
    for(i in 1:length(maxl)){
      if(length(arg) > maxl[i]){
        nameMax = maxl[i]
        if(!is.null(names(maxl[i])))
          nameMax = names(maxl[i])
        stop(paste("argument", substitute(arg) ,"is greater then", nameMax))
      }
    }
} 