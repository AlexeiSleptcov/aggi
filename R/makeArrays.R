#' Make arrays
#'
#' Reassambles arrays from RGList
#' 
#' @param data RGList only
#' 
#' @export

makeArrays <- function(data){
  if(!inherits(data, c("RGList", "EListRaw")))
    stop("Data must RGList class object",  call. = FALSE)
  
  # find missing values
  temp   = data.frame(row=data$genes$Row, col=data$genes$Col, stringsAsFactors = FALSE)
  delete = which(is.na(melt(acast(temp, col~row, value.var="col"))$value))
  
  if(inherits(data,"RGList")){
    # row, column, array, 1=R/2=G
    result = array(NA, c(max(data$genes$Row), max(data$genes$Col), ncol(data), 2))
    for(n in 1:ncol(data)){
      gridR = data.frame(row=data$genes$Row, col=data$genes$Col,
                         int=data$R[,n], stringsAsFactors = FALSE)
      gridG = data.frame(row=data$genes$Row, col=data$genes$Col,
                         int=data$G[,n], stringsAsFactors = FALSE)
      result[,,n,1] = acast(gridR, row~col, value.var="int")
      result[,,n,2] = acast(gridG, row~col, value.var="int")
    }
  } else {
    result = array(NA, c(max(data$genes$Row), max(data$genes$Col), ncol(data)))
    # row, column, array
    for(n in 1:ncol(data)){
      gridE = data.frame(row=data$genes$Row, col=data$genes$Col,
                         int=data$E[,n], stringsAsFactors = FALSE)
      result[,,n] = acast(gridE, row~col, value.var="int")
    }
  }
  
 # descriptives
  
  returns = list()
  returns[["Arrays"]] = result
  controls = data.frame(row=data$genes$Row, col=data$genes$Col,
                        int=data$genes$ControlType, stringsAsFactors = FALSE)
  returns[["ControlType"]] = acast(controls, row~col, value.var="int")
  returns[["MissingString"]] = delete
  
  if(inherits(data,"RGList")){
    class(returns) = c("RGarray", "list")
  } else {
    class(returns) = c("Earray", "list")
  }
  returns
}