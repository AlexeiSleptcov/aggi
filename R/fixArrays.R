#' Slope correction intra-array
#' 
#' Correction slope of arrays by Negatives probe and based on Affine transformation. 
#' Optimal span for 2DLoess normalization by used generalized cross-validation. See 
#' \code{gcv} in \link{\code{fANCOVA}}.
#'  
#' @param data RGList, EListRaw
#' 
#' @export

fixArrays <- function(data){
  
  ## find optimal span
  # optimal span by fANCOVA
  # bias-corrected Akaike information criterion (aicc)
  # and generalized cross-validation (gcv)
  ##
  optimal.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
    as.crit <- function(x) {
      span   <- x$pars$span
      traceL <- x$trace.hat
      sigma2 <- sum(x$residuals^2)/(x$n - 1)
      aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2)
      gcv    <- x$n * sigma2/(x$n - traceL)^2
      result <- list(span = span, aicc = aicc, gcv = gcv)
      return(result)
    }
    criterion <- match.arg(criterion)
    fn <- function(span) {
      mod <- update(model, span = span)
      as.crit(mod)[[criterion]]
    }
    result <- optimize(fn, span.range)
    return(list(span = result$minimum, criterion = result$objective))
  }
  
  # Greatest Common Divisor
  gcd <- function(x,y) {
    r <- x%%y;
    return(ifelse(r, gcd(y, r), y))
  }
  
  # RUN
   if(inherits(data, "RGList")){
     for(smp in 1:ncol(data)){
       
       dneg = data
       adata = makeArrays(data)
       
       for(cyan in c("R", "G")){
         if(cyan == "R"){
           neg = x = adata$Arrays[,,smp, 1]
         } else {
           neg = x = adata$Arrays[,,smp, 2]
         }
         dneg[[cyan]][data$genes$ControlType != -1,] = NA
         neg[data$ControlType != -1] = NA
         
         # smoothing NegativesE
         
         sneg = makeSmall(neg)
         
         dneg[[cyan]][which(data$genes$Row == 1 & data$genes$Col == 1)] <- 
           median(sneg[1:5,1:5], na.rm = TRUE) # upRow leftCol
         dneg[[cyan]][which(data$genes$Row == 1 & data$genes$Col == max(data$genes$Col))] <- 
           median(sneg[(dim(sneg)[1]-5):dim(sneg)[1],1:5], na.rm = TRUE) # upRow rightCol
         dneg[[cyan]][which(data$genes$Row == max(data$genes$Row) & data$genes$Col == max(data$genes$Col))] <- 
           median(sneg[(dim(sneg)[1]-5):dim(sneg)[1],(dim(sneg)[2]-5):dim(sneg)[2]], na.rm = TRUE) # downRow rightCol
         dneg[[cyan]][which(data$genes$Row == max(data$genes$Row) & data$genes$Col == 1)] <- 
           median(sneg[1:5,(dim(sneg)[2]-5):dim(sneg)[2]], na.rm = TRUE) # downRow leftCol
         
         
         moddf = data.frame(Int = dneg[[cyan]][,smp], Row = data$genes$Row, Col = data$genes$Col, 
                            stringsAsFactors = FALSE)
         model = loess(Int ~ Row * Col, 
                       degree = 1, span = 1, family = "symmetric", data = moddf)
         moddf$Fit = predict(model, moddf[,-1])
         
         data[[cyan]][,smp] = data[[cyan]][,smp] - moddf$Fit
         
         # smoothing Regular
         sx = makeSmall(x)
         sxm = na.omit(melt(sx)[,-4])
         names(sxm) = c("Row", "Col", "Int")
         
         xgcd = gcd(dim(x)[1], dim(x)[2])
         sxm$Row = c(1,(sxm$Row * xgcd)[-1])
         sxm$Col = c(1,(sxm$Col * xgcd)[-1])
         
         optmodel = loess(Int ~ Row * Col, 
                          degree = 1, span = 0.1, family = "symmetric", data = sxm)
         optSpan <- optimal.span(optmodel, criterion = "gcv")$span
         model = loess(Int ~ Row * Col, 
                       degree = 1, span = optSpan, family = "symmetric", data = sxm)
         
         slope = predict(model, data.frame(Row = data$genes$Row, Col = data$genes$Col))
         
         d = data[[cyan]][,smp]*(median(slope, na.rm = TRUE)/slope)
         data[[cyan]][,smp] = d - summary(d)[1] + 25 
       }
     }
     
   } else if (inherits(data, "EListRaw")){
     
     dneg = data
     dneg$E[data$genes$ControlType != -1,] = NA
     
     adata = makeArrays(data)
     
     for(smp in 1:ncol(data)){
       neg = x = adata$Arrays[,,smp]
       neg[data$ControlType != -1] = NA
       
       # smoothing Negatives
       
       sneg = makeSmall(neg)
       
       dneg$E[which(data$genes$Row == 1 & data$genes$Col == 1)] <- 
         median(sneg[1:5,1:5], na.rm = TRUE) # upRow leftCol
       dneg$E[which(data$genes$Row == 1 & data$genes$Col == max(data$genes$Col))] <- 
         median(sneg[(dim(sneg)[1]-5):dim(sneg)[1],1:5], na.rm = TRUE) # upRow rightCol
       dneg$E[which(data$genes$Row == max(data$genes$Row) & data$genes$Col == max(data$genes$Col))] <- 
         median(sneg[(dim(sneg)[1]-5):dim(sneg)[1],(dim(sneg)[2]-5):dim(sneg)[2]], na.rm = TRUE) # downRow rightCol
       dneg$E[which(data$genes$Row == max(data$genes$Row) & data$genes$Col == 1)] <- 
         median(sneg[1:5,(dim(sneg)[2]-5):dim(sneg)[2]], na.rm = TRUE) # downRow leftCol
       
       
       moddf = data.frame(Int = dneg$E[,smp], Row = data$genes$Row, Col = data$genes$Col, 
                          stringsAsFactors = FALSE)
       model = loess(Int ~ Row * Col, 
                     degree = 1, span = 1, family = "symmetric", data = moddf)
       moddf$Fit = predict(model, moddf[,-1])
       
       data$E[,smp] = data$E[,smp] - moddf$Fit
       
       # smoothing Regular
       sx = makeSmall(x)
       sxm = na.omit(melt(sx)[,-4])
       names(sxm) = c("Row", "Col", "Int")
       
       xgcd = gcd(dim(x)[1], dim(x)[2])
       sxm$Row = c(1,(sxm$Row * xgcd)[-1])
       sxm$Col = c(1,(sxm$Col * xgcd)[-1])
       
       optmodel = loess(Int ~ Row * Col, 
                        degree = 1, span = 0.1, family = "symmetric", data = sxm)
       optSpan <- optimal.span(optmodel, criterion = "gcv")$span
       model = loess(Int ~ Row * Col, 
                     degree = 1, span = optSpan, family = "symmetric", data = sxm)
       
       slope = predict(model, data.frame(Row = data$genes$Row, Col = data$genes$Col))
       
       d = data$E[,smp]*(median(slope, na.rm = TRUE)/slope)
       data$E[,smp] = d - summary(d)[1] + median(moddf$Fit, na.rm = TRUE) 
     }
     
   } else stop("Data must be RGList or EListRaw", call. = FALSE)
  data
}


