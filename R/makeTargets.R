#' Make targets object
#' 
#' Create targets object based on TXT files localized in filepath.
#' 
#' @param path filepath, where content only your TXT FE-files. Default: \code{getwd()}
#' @param manames names for your arrays. If argument is \code{NULL} (default) names will
#' be assigned a serial number, for instance: MA1, MA2 etc.
#' 
#' @example mytargets <- makeTargets()
#' 
#' @export

makeTargets <- function(path = getwd(), 
                        manames = NULL){
  
  files = list.files(filepath, "*.txt")
  
  if(is.null(manames)){
    manames = paste("MA", 1:length(files), sep = "")
  }
  
  if(length(files) != length(manames))
    stop("Length files and manames are different", call. = FALSE)
  
  data.frame(Name = manames, FileName = files, stringsAsFactors = FALSE)
}