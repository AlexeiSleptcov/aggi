#' Noise reduction for aCGH
#'
#' Noise reduction only for array CGH data
#' 
#' @param data RGList only
#' @param N numeric, i
#' @param method iqr or mad
#' 
#' @export

noiseReductionCGH <- function(data, N = 200, method = c("mad", "iqr")){
  if(!inherits(data, "RGList"))
    stop("First, data must be RGList", call. = FALSE)
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
    data$R[,i] = ( ( (Rf2-Rm)*(N/Rmad) )+Rm )
    data$G[,i] = ( ( (Gf2-Gm)*(N/Gmad) )+Gm )
    data$R[,i] = data$R[,i]-median(data$R[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)
    data$G[,i] = data$G[,i]-median(data$G[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)
    
    
    cmeanNeg = ((median(data$R[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)+median(data$G[data$genes$SubTypeMask == 1028,i], na.rm = TRUE)))/2
    cmean = ((median(data$R[,i], na.rm = TRUE)+median(data$G[,i], na.rm = TRUE)))/2
    
    Rz = data$R[,i]-cmeanNeg
    Gz = data$G[,i]-cmeanNeg
    
    data$R[,i] = (data$R[,i]*(cmean/median(data$R[,i], na.rm = TRUE)))+cmeanNeg
    data$G[,i] = (data$G[,i]*(cmean/median(data$G[,i], na.rm = TRUE)))+cmeanNeg
  }
  data
}
