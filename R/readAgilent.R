#' Read agilent files
#'
#' Read Feature Extraction's txt files.
#'
#' @param targets Feature Extraction's files. TXT format files only.
#' @param path filepath, where FE-files are stored. Default: getwd().
#' @param format which format your FE-files have? Default: \code{"compact"}
#' @param macolors read two-color (2) or one-color (1) microarray.
#' @param collect what the data signals you want to process? Default: \code{"processed"}.
#'
#' @return Read Feature Extraction's files.
#' @return Collecting data from the format \code{"compact"} can be taken on the processed signals
#' (\code{"processed"}), raw data (\code{"mean"} or \code{"median"}) or LogRatio (\code{"LogRatio"}).
#' At the same time, if format selected as \code{"full"} you may collect normalized data (\code{"DyeNorm"}).
#'
#' @examples 
#' # targets <- read.table("targets.txt", header = TRUE, stringsAsFactors = FALSE)
#' # MyArrays = readAgilent(targets, collect = "median")
#'
#' @export

readAgilent <- function(targets = NULL,
                        path = getwd(),
                        format = c("compact", "full"),
                        macolors = 2,
                        collect = c("processed", "median", "mean", "LogRatio", "DyeNorm")
                        ){
  # check paramethers
  if(is.null(targets))
    stop("First of all, needs targets object", call. = FALSE)

  collect = match.arg(collect)
  format = match.arg(format)

  if(format == "compact" & collect == "DyeNorm")
    stop("Argument DyeNorm only in full format" , call. = FALSE)

  if(macolors == 1){
    if(collect == "LogRatio")
      stop("Argument LogRatio only in two-color microarray" , call. = FALSE)

    columns = c(E="gMedianSignal", Emean="gMeanSignal")
    switch(collect,
           "processed" = {
             columns = list(E="gProcessedSignal")
           },
           "median" = {
             columns = list(E="gMedianSignal", Emean="gMeanSignal")
           },
           "mean" = {
             columns = list(E="gMeanSignal", Emean="gMedianSignal")
           },
           "DyeNorm" = {
             columns = list(E="gDyeNormSignal")
           })
    other.columns = list(gIsFeatPopnOL="gIsFeatPopnOL",
                         gIsSaturated="gIsSaturated",
                         gIsFeatNonUnifOL="gIsFeatNonUnifOL")

  } else if( macolors == 2) {
    # constituitive paramethers
    other = list(flag = "flag",
                 annotation = "annotation",
                 gIsFeatNonUnifOL="gIsFeatNonUnifOL",
                 rIsFeatNonUnifOL="rIsFeatNonUnifOL",
                 gIsBGNonUnifOL="gIsBGNonUnifOL",
                 rIsBGNonUnifOL="rIsBGNonUnifOL",
                 gIsFeatPopnOL="gIsFeatPopnOL",
                 rIsFeatPopnOL="rIsFeatPopnOL",
                 gIsBGPopnOL="gIsBGPopnOL",
                 rIsBGPopnOL="rIsBGPopnOL",
                 rIsSaturated="rIsSaturated",
                 gIsSaturated="gIsSaturated")

    # Description
    # FeatNonUnifOL = signal is non-uniformity outlier
    # FeatPopnOL = signal is population outlier
    # PosAndSignif = signal exceeds significantly background
    # WellAboveBG = signal exceeds more significantly background

    switch(collect,
           "DyeNorm" = {
             columns = list(R="rDyeNormSignal",G="gDyeNormSignal")
             other.columns = c(other, "rMeanSignal", "gMeanSignal",
                               "rMedianSignal", "gMedianSignal")
           },
           "ProcessedSignal" = {
             columns = list(R="rProcessedSignal",G="gProcessedSignal",
                            Rb = "rProcessedBackground", Gb = "gProcessedBackground")
             other.columns = c(other, "rMeanSignal", "gMeanSignal",
                               "rMedianSignal", "gMedianSignal")
           },
           "median" = {
             columns = list(R = "rMedianSignal", G = "gMedianSignal",
                            Rb = "rBGMedianSignal", Gb = "gBGMedianSignal")
             other.columns = c(other, "rMeanSignal", "gMeanSignal")
           },
           "mean" = {
             columns = list(R = "rMeanSignal", G = "gMeanSignal",
                            Rb = "rBGMeanSignal", Gb = "gBGMeanSignal")
             other.columns = c(other, "rMedianSignal", "gMedianSignal")
           },
           "LogRatio" = {
             logs = list(LR="LogRatio",LRE="LogRatioError")
             other.columns = c(logs, other)
             columns = NULL
           })
  } else stop("Argument macolors might be only 1 or 2", call. = FALSE)

  read.data = read.maimages(files,
                            source = "agilent.median",
                            path=path,
                            names = basename(targets$Name),
                            columns=columns,
                            other.columns=other.columns,
                            annotation = c("Row", "Col", "Sequence", "ControlType",
                                           "ProbeName", "GeneName", "SystematicName",
                                           "identifier", "SubTypeMask"))
  read.data
}
