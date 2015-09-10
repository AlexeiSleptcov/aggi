#' Find outliners in arrays
#'
#' Exclude outliners by 2 ways, either by FE outliner flags or find spost with unequal 
#' ratio of mean and median pixel intensity.
#' 
#' @param data RGList, EListRaw
#' @param smp numeric, number of sample
#' @param flags logical, if \code{TRUE} find outliners by flags else by ratio of mean and median
#' @param verbose logical, if \code{TRUE} create plot of delta (only if delta is \code{TRUE})
#' 
#' @export

findOutliners <- function(data, smp = 1, flags = TRUE, verbose = FALSE){
  
  if(!inherits(data,c("RGList", "EListRaw")))
    stop("First, data must be RGList or EListRaw", call.=FALSE)
  
  x <- data$other
  results = list()
  if(flags){
    if(inherits(data,"RGList")){
      w = list()
      # FeatNonUnifOL = signal is non-uniformity outlier
      w[[1]] = which(x$gIsFeatNonUnifOL[,smp] == 1 | x$rIsFeatNonUnifOL[,smp] == 1)
      # FeatPopnOL = signal is population outlier
      w[[2]] = which(x$gIsFeatPopnOL[,smp]    == 1 | x$rIsFeatPopnOL[,smp]    == 1)
      # PosAndSignif = signal exceeds significantly background
      w[[3]] = which(x$gIsSaturated[,smp]     == 1 | x$rIsSaturated[,smp]     == 1)
      # WellAboveBG = signal exceeds more significantly background
      probes = unlist(w)
      probes[!duplicated(probes)]
    } else {
      w = list()
      # FeatNonUnifOL = signal is non-uniformity outlier
      w[[1]] = which(x$gIsFeatNonUnifOL[,smp] == 1)
      # FeatPopnOL = signal is population outlier
      w[[2]] = which(x$gIsFeatPopnOL[,smp]    == 1)
      # PosAndSignif = signal exceeds significantly background
      w[[3]] = which(x$gIsSaturated[,smp]     == 1)
      # WellAboveBG = signal exceeds more significantly background
      probes = unlist(w)
      probes[!duplicated(probes)]
    }
    
  } else {
    # remove "bad spots" = spost with unequal ratio of mean and
    #  median pixel intensity
    bs = c()
    if(inherits(data,"RGList")){
      for(li in c("R","G")){
        if (!is.null(x$rMedianSignal) && !is.null(x$rMeanSignal)){
          if(li == "R"){
            X  = x$rMedianSignal[,smp]
            Xm = x$rMeanSignal[,smp]
          } else {
            X  = x$gMedianSignal[,smp]
            Xm = x$gMeanSignal[,smp]
          }
        } else {
          X  = data[[li]][,smp]
          if(li == "R"){
            if(is.null(x$rMedianSignal)){
              Xm = x$rMeanSignal[,smp]
            } else {
              Xm = x$rMedianSignal[,smp]
            } # Rm signal
          } else {
            if(is.null(x$rMedianSignal)){
              Xm = x$gMeanSignal[,smp]
            } else {
              Xm = x$gMedianSignal[,smp]
            } # Gm signal
          }
        } # R or G signal
        
        d = log2(X)-log2(Xm)
        s = log2(X)
        # the relevant log ratio depends on the signal intensity
        # very low intensity spots will likely have larger deltas
        
        if(verbose){
          plot(s, d, pch=".", main=paste("Bad Spots in", colnames(data$R)[smp]), 
               xlab=paste("Log2", li, "signal"), ylab="Delta")
          # draw the cut-of line as an exponentila function of the signal
          x = seq(0,16,len=100)
          y = 8*2^(-x/2) + 0.25
          lines(x,y,col="red"); lines(x,-y,col="red")
          # identify the spots with too large deltas, compared to other
          #  spots of similar signal intensity
          bad = abs(d) > (8*2^(-s/2) + 0.25)
          # mark theses points in the plot
          points(s[bad], d[bad], pch=16, col="red")
        } else {
          bad = abs(d) > (8*2^(-s/2) + 0.25)
        }
        bs = c(bs, which(bad))
      }
    } else {
      X  = data$E[,smp]
      if(is.null(data$Emean))
        stop("Emean is NULL", call. = FALSE)
      Xm = data$Emean[,smp]
      
      d = log2(X)-log2(Xm)
      s = log2(X)
      # the relevant log ratio depends on the signal intensity
      # very low intensity spots will likely have larger deltas
      
      if(verbose){
        plot(s, d, pch=".", main=paste("Bad Spots in", colnames(data$E)[smp]), 
             xlab=paste("Log2 signals"), ylab="Delta")
        # draw the cut-of line as an exponentila function of the signal
        x = seq(0,16,len=100)
        y = 8*2^(-x/2) + 0.25
        lines(x,y,col="red"); lines(x,-y,col="red")
        # identify the spots with too large deltas, compared to other
        #  spots of similar signal intensity
        bad = abs(d) > (8*2^(-s/2) + 0.25)
        # mark theses points in the plot
        points(s[bad], d[bad], pch=16, col="red")
      } else {
        bad = abs(d) > (8*2^(-s/2) + 0.25)
      }
      bs = c(bs, which(bad))
    }
    bs[!duplicated(bs)]
  }
}