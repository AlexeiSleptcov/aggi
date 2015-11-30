#' Multicore Circular Binary Segmentation
#' 
#' Detect outliers and smooth the data then perform circular binary segmentation.
#' 
#' @param data MAList
#' @param sampleid sample identifier. If missing the samples are named by prefixing "Sample" to consecutive integers.
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
#' @return See details \code{\link{CNA}}, \code{\link{smooth.CNA}}, \code{\link{segmByCluster}}, \code{\link{segment}}
#' 
#' @export


runCBS <- function(data, 
                   sampleid = NULL,
                   weights = NULL, 
                   alpha = 0.01, 
                   nperm = 10000, 
                   p.method = c("hybrid", "perm"), 
                   min.width = 3,
                   kmax = 25, 
                   nmin = 200, 
                   eta = 0.05,
                   sbdry = NULL, 
                   trim = 0.025, 
                   undo.splits = c("sdundo", "prune", "none"), 
                   undo.prune = 0.05, 
                   undo.SD = 3, 
                   verbose = 1,
                   cluster=NULL){
  
  if(!inherits(data, "MAList"))
    stop("Data must be MAList", call. = FALSE)
  
  #cat("GC correction...\n") 
  #corr <- GCcorrection.agi(data$M, cluster=cluster)
  genomdat = chrmap(data, obj = FALSE, verbose = 0, ingenes = TRUE)
  if(is.null(sampleid)) sampleid = colnames(data)
  a = CNA(genomdat$M, genomdat$genes$chr, genomdat$genes$maploc, sampleid = sampleid)
  #rm(list=c("maploc","chrom"))
  x = smooth.CNA(a, smooth.region=10)
  if(is.null(cluster)){
    segment(x,alpha = alpha, 
                  nperm = nperm, 
                  p.method = match.arg(p.method), 
                  min.width = min.width,
                  kmax = kmax, 
                  nmin = nmin, 
                  eta = eta,
                  sbdry = sbdry, 
                  trim = trim, 
                  undo.splits = match.arg(undo.splits), 
                  undo.prune = undo.prune, 
                  undo.SD = undo.SD, 
                  verbose = verbose)
  } else {
    segmByCluster(x, weights = weights, 
                  alpha = alpha, 
                  nperm = nperm, 
                  p.method = match.arg(p.method), 
                  min.width = min.width,
                  kmax = kmax, 
                  nmin = nmin, 
                  eta = eta,
                  sbdry = sbdry, 
                  trim = trim, 
                  undo.splits = match.arg(undo.splits), 
                  undo.prune = undo.prune, 
                  undo.SD = undo.SD, 
                  verbose = verbose,
                  cluster = cluster)
  }
   
  
}
