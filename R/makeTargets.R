#' Make targets object
#' 
#' Create targets object based on files localized in filepath.
#' 
#' @param path filepath, where content only your files. Default: \code{getwd()}
#' @param manames names for your arrays. If argument is \code{NULL} (default) names will
#' @param fileformat Default \code{"*.txt"}
#' be assigned a serial number, for instance: MA1, MA2 etc.
#' 
#' @examples # mytargets <- makeTargets()
#' 
#' @export

makeTargets <- function(path = getwd(), 
                        manames = NULL,
                        fileformat = "*.txt"){
  
  files = list.files(path, fileformat)
  
  if(is.null(manames)){
    manames = paste("MA", 1:length(files), sep = "")
  }
  
  if(length(files) != length(manames))
    stop("Length files and manames are different", call. = FALSE)
  
  data.frame(names = manames, filenames = files, stringsAsFactors = FALSE)
}