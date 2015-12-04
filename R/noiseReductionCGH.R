#' Noise reduction for aCGH
#'
#' Noise reduction only for array CGH data
#' 
#' @param data RGList only
#' @param N numeric, i
#' @param method iqr or mad
#' 
#' @export

noiseReductionCGH <- function(data, N = NULL, method = c("mad", "iqr")){
  if(!inherits(data, c("RGList", "list")))
    stop("First, data must be RGList", call. = FALSE)
  if(inherits(data, "list")){
    if(!(exists("R", where = data) && exists("G", where = data)))
      stop("Data do not have R and G objects", call. = FALSE)
  }
  
  for(i in 1:ncol(data$R)){
    Rf2 = data$R[,i]
    Gf2 = data$G[,i]
    Rm = median(Rf2, na.rm = TRUE)
    Gm = median(Gf2, na.rm = TRUE)
    if(match.arg(method) == "mad") {
      Rmad = mad(Rf2, na.rm = TRUE)/2
      Gmad = mad(Gf2, na.rm = TRUE)/2
    } else {
      Rmad = IQR(Rf2, na.rm = TRUE)/2
      Gmad = IQR(Gf2, na.rm = TRUE)/2
    }
    if(is.null(N))
      N = mean(c(Rmad, Gmad))
    
    # AFFINE TRANSFORMATION
    # SUBTRACT MEDIAN BEFORE & ADD AFTER
    # FITTING MAD TO COMMON MAD
    data$R[,i] = ( ( (Rf2-Rm)*(N/Rmad) )+Rm )
    data$G[,i] = ( ( (Gf2-Gm)*(N/Gmad) )+Gm )
    
#     if(inherits(data, "RGlist")){
#       data$R[,i] = data$R[,i]-median(data$R[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)
#       data$G[,i] = data$G[,i]-median(data$G[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)
#     }
    
    # if(inherits(data, "list")){
    
    # FIT COMMON MEDIAN
      cmean = ((median(data$R[,i], na.rm = TRUE)+median(data$G[,i], na.rm = TRUE)))/2
      
      data$R[,i] = data$R[,i]*(cmean/median(data$R[,i], na.rm = TRUE))
      data$G[,i] = data$G[,i]*(cmean/median(data$G[,i], na.rm = TRUE))
#     } else {
#       cmeanNeg = ((median(data$R[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)+median(data$G[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)))/2
#       cmean = ((median(data$R[,i], na.rm = TRUE)+median(data$G[,i], na.rm = TRUE)))/2
#       
#       Rz = data$R[,i]-cmeanNeg
#       Gz = data$G[,i]-cmeanNeg
#       
#       data$R[,i] = (data$R[,i]*(cmean/median(data$R[,i], na.rm = TRUE)))+cmeanNeg
#       data$G[,i] = (data$G[,i]*(cmean/median(data$G[,i], na.rm = TRUE)))+cmeanNeg
#     }
  }
  data
}

