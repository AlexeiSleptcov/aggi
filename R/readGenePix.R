#' Read genepix files
#'
#' Read GenePix grp files.
#'
#' @param targets GenePix's files. GPR format files only.
#' @param path filepath, where GP-files are stored. Default: getwd().
#' @param custom.data custom data. For example \code{c(E = "F532 Median", Eb = "F532 MedianBG")}
#' @param custom.ann custom annotation. For example \code{c(Row="Row", Col="Column")}
#'
#' @details Targets object must be class \code{data.frame}, and have two required column \code{name} 
#' (eg.: "array1", "array2") and \code{filenames} (eg. "2802860_12158125.gpr").
#'  
#' @seealso \code{\link{readAgilent}}  
#'  
#' @examples 
#' # targets <- read.table("targets.txt", header = TRUE, stringsAsFactors = FALSE)
#' # MyArrays <- readGenePix(targets)
#'
#' @export

readGenePix <- function(targets = NULL,
                        path = getwd(),
                        custom.data = NULL,
                        custom.ann = NULL){
  # check paramethers
  if(is.null(targets))
    stop("First of all, needs targets object", call. = FALSE)
  
  if(format == "compact" & collect == "DyeNorm")
    stop("Argument DyeNorm only in full format" , call. = FALSE)
  
  columns.data = c(E="F532 Median",
                   Emean="F532 Mean")
  
  columns.annotation = c(Row="Row", 
                         Col="Column",
                         ProbeName="ID", 
                         ControlType="ControlType",
                         GeneName="GeneName")
  
  read.data = read.maimages(targets$filenames,
                            source = "genepix.median",
                            path=path,
                            names = basename(targets$names),
                            columns=columns,
                            other.columns=other.columns,
                            annotation = columns.annotation)
  read.data
}
