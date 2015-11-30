#' Multicore Genome Segmentation Program 
#' 
#' This program segments DNA copy number 
#' data into regions of estimated 
#' equal copy number using circular binary segmentation (CBS) (\code{\link{DNAcopy}} package)
#' and multicore (\code{\link{doParallel}}).
#' 
#' @param x	an object of class CNA
#' @param weights	a vector of weights for the probes. The weights should be inversely proportional to their variances. Currently all weights should be positive i.e. remove probes with zero weight prior to segmentation.
#' @param alpha	significance levels for the test to accept change-points.
#' @param nperm	number of permutations used for p-value computation.
#' @param p.method	method used for p-value computation. For the "perm" method the p-value is based on full permutation. For the "hybrid" method the maximum over the entire region is split into maximum of max over small segments and max over the rest. Approximation is used for the larger segment max. Default is hybrid.
#' @param min.width	 the minimum number of markers for a changed segment. The default is 2 but can be made larger. Maximum possible value is set at 5 since arbitrary widths can have the undesirable effect of incorrect change-points when a true signal of narrow widths exists.
#' @param kmax	 the maximum width of smaller segment for permutation in the hybrid method.
#' @param nmin	 the minimum length of data for which the approximation of maximum statistic is used under the hybrid method. should be larger than 4*kmax
#' @param eta	 the probability to declare a change conditioned on the permuted statistic exceeding the observed statistic exactly j (= 1,...,nperm*alpha) times.
#' @param sbdry	 the sequential boundary used to stop and declare a change. This boundary is a function of nperm, alpha and eta. It can be obtained using the function "getbdry" and used instead of having the "segment" function compute it every time it is called.
#' @param trim	 proportion of data to be trimmed for variance calculation for smoothing outliers and undoing splits based on SD.
#' @param undo.splits	 A character string specifying how change-points are to be undone, if at all. Default is "none". Other choices are "prune", which uses a sum of squares criterion, and "sdundo", which undoes splits that are not at least this many SDs apart.
#' @param undo.prune	 the proportional increase in sum of squares allowed when eliminating splits if undo.splits="prune".
#' @param undo.SD	 the number of SDs between means to keep a split if undo.splits="sdundo".
#' @param verbose	 level of verbosity for monitoring the program's progress where 0 produces no printout, 1 prints the current sample, 2 the current chromosome and 3 the current segment. The default level is 1.
#' @param cluster A specification appropriate to the type of cluster.
#' 
#' @return See details \code{\link{segment}} and \code{\link{makeCluster}}
#' 
#' @export

segmByCluster <- function (x, 
                           weights = NULL, 
                           alpha = 0.01, 
                           nperm = 10000, 
                           p.method = c("hybrid", "perm"), 
                           min.width = 2,
                           kmax = 25, 
                           nmin = 200, 
                           eta = 0.05,
                           sbdry = NULL, 
                           trim = 0.025, 
                           undo.splits = c("none", "prune", "sdundo"), 
                           undo.prune = 0.05, 
                           undo.SD = 2, 
                           verbose = 1,
                           cluster=getDoParWorkers()){
  if (!inherits(x, "CNA")) 
    stop("First arg must be a copy number array object")
  call <- match.call()
  if (min.width < 2 | min.width > 5) 
    stop("minimum segment width should be between 2 and 5")
  if (nmin < 4 * kmax) 
    stop("nmin should be >= 4*kmax")
  if (missing(sbdry)) {
    if (nperm == 10000 & alpha == 0.01 & eta == 0.05) {
      if (!exists("default.DNAcopy.bdry")) 
        data(default.DNAcopy.bdry, package = "DNAcopy", 
             envir = environment())
      sbdry <- default.DNAcopy.bdry
    }
    else {
      max.ones <- floor(nperm * alpha) + 1
      sbdry <- getbdry(eta, nperm, max.ones)
    }
  }
  weighted = !is.null(weights)
  if (weighted) {
    if (length(weights) != nrow(x)) 
      stop("length of weights should be the same as the number of probes")
    if (min(weights) <= 0) 
      stop("all weights should be positive")
  }
  sbn <- length(sbdry)
  nsample <- ncol(x) - 2
  sampleid <- colnames(x)[-(1:2)]
  uchrom <- unique(x$chrom)
  data.type <- attr(x, "data.type")
  p.method <- match.arg(p.method)
  undo.splits <- match.arg(undo.splits)
  segres <- list(data = x)
  allsegs <- list(ID = NULL, chrom = NULL, loc.start = NULL, loc.end = NULL, num.mark = NULL,seg.mean = NULL)
  segRows <- list(startRow = NULL, endRow = NULL)
  if (verbose >= 1) {
    pb <- txtProgressBar(min = 0, max = nsample, style = 3)
  }
  for (isamp in 1:nsample) {
    genomdati <- x[, isamp + 2]
    ina <- which(is.finite(genomdati))
    genomdati <- genomdati[ina]
    trimmed.SD <- sqrt(DNAcopy:::trimmed.variance(genomdati, trim))
    chromi <- x$chrom[ina]
    if (weighted) {
      wghts <- weights[ina]
    }
    else {
      wghts <- NULL
    }
    sample.lsegs <- NULL
    sample.segmeans <- NULL
    
    if(require(doParallel)){ cl <- makeCluster(cluster); registerDoParallel(cl)}
    segci <<- foreach(ic=uchrom, .combine=rbind) %dopar% DNAcopy:::changepoints(genomdati[chromi == ic], data.type, 
                                                                                alpha, wghts, sbdry, sbn, nperm, p.method, min.width, kmax, nmin, trimmed.SD, undo.splits, undo.prune, undo.SD, verbose)
    if(require(doParallel)) stopCluster(cl)
    sample.lsegs = unname(unlist(segci[,1]))
    sample.segmeans = unname(unlist(segci[,2]))
    
    sample.nseg <- length(sample.lsegs)
    sample.segs.start <- ina[cumsum(c(1, sample.lsegs[-sample.nseg]))]
    sample.segs.end <- ina[cumsum(sample.lsegs)]
    allsegs$ID <- c(allsegs$ID, rep(isamp, sample.nseg))
    allsegs$chrom <- c(allsegs$chrom, x$chrom[sample.segs.end])
    allsegs$loc.start <- c(allsegs$loc.start, x$maploc[sample.segs.start])
    allsegs$loc.end <- c(allsegs$loc.end, x$maploc[sample.segs.end])
    allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
    allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
    segRows$startRow <- c(segRows$startRow, sample.segs.start)
    segRows$endRow <- c(segRows$endRow, sample.segs.end)
    if (verbose >= 1) {
      setTxtProgressBar(pb, isamp)
    }
  }
  
  if (verbose >= 1) {
    close(pb)
  }
  allsegs$ID <- sampleid[allsegs$ID]
  allsegs$seg.mean <- round(allsegs$seg.mean, 4)
  allsegs <- as.data.frame(allsegs)
  allsegs$ID <- as.character(allsegs$ID)
  segres$output <- allsegs
  segres$segRows <- as.data.frame(segRows)
  segres$call <- call
  if (weighted) {
    segres$weights <- weights
  }
  class(segres) <- "DNAcopy"
  closeAllConnections()
  segres
}