#' Make small array 
#'
#' The reduction is based on the Greatest Common Divisor and median
#' 
#' @param data matrix
#' 
#' @return The size of reduced matrix is determined by the Greatest Common Divisor. 
#' Reduced values will be aligned by median.
#'
#' @keywords internal
#' @export


makeSmall <- function(data){
  
  # http://stackoverflow.com/questions/21502181/finding-the-gcd-without-looping-r
  gcd <- function(x,y) {
    r <- x%%y;
    return(ifelse(r, gcd(y, r), y))
  }
  
  square = gcd(dim(data)[1], dim(data)[2])
  
  n.row = nrow(data)
  n.col = ncol(data)
  
  s.row = 1:round(n.row/square)
  s.col = 1:round(n.col/square)
  
  d.row = 1:n.row
  d.col = 1:n.col
  
  rem.r = n.row %/% square
  rem.c = n.col %/% square
  
  s.grid = matrix(NA, ncol= max(s.col), nrow = max(s.row))
  
  pr = 1
  for(ir in s.row){
    if(ir > rem.r){
      itr = rem.r
    } else {
      itr = ir
    }
    
    l.row = d.row[pr]:(d.row[(itr*square)])
    pc = 1
    for(ic in s.col){
      if(ic > rem.c) {
        itc = rem.c
      } else {
        itc = ic
      }
      
      l.col = d.col[pc]:(d.col[(itc*square)]) 
      
      s.grid[ir,ic] = median(data[l.row, l.col], na.rm=TRUE)
      
      pc = ic*square+1
    } # end for column
    pr = ir*square+1
  } # end for row
  s.grid
}