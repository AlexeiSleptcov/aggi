#' Plot SNP distribution
#'
#' @param data MAList or RGList with SNP probes
#' @param array number of array
#' @param reference matrix, data of cut/uncut alleles
#' @param type type of SNP ditribution: \code{log} or \code{theta}
#' @param SNPprobe numeric, type number of SNP probes. Default: -15000
#' 
#' @details Argument \code{array} might be numeric, sequence number of array or 
#' parametre \code{"all"} to plot all samples in one plot.
#' @details If number of \code{SNPprobe} not found in \code{ControlType} in arrays
#' plot choose all probes.
#' 
#' @export

gplotSNP <- function(data, 
                     array = 1, 
                     reference = c("none", "male", "female"), 
                     type = c("log","theta"),
                     SNPprobe = -15000){
  
  if(inherits(data,"MAList")){
  } else if (inherits(data,"RGList")){
    data = MA.RG(data)
  } else stop("Data must be MAList or RGList", call.=FALSE)
  
  if( any(unique(data$genes$ControlType) %in% SNPprobe) ){
    data = data[data$genes$ControlType == SNPprobe,]
  } else warning("SNP probes not found", call. = FALSE)
   
  if(array == "all"){
    dataSNP = data.frame(ID = 1:nrow(data$M), data$M, stringsAsFactors = FALSE)
    dataSNP = melt(dataSNP, na.rm = TRUE, id = "ID")
    names(dataSNP) = c("ID", "samples", "M")
    mediana = median(dataSNP$M, na.rm = TRUE)
    mytitle = "Distribution of SNP probes"
  } else { 
    mediana = median(data$M[,array], na.rm=TRUE)
    mytitle = paste("Distribution of SNP probes in", colnames(data$M)[array])
  }
  
  reference = match.arg(reference)
  if(reference != "none"){
    # devtools::use_data(AgilentEuroMale, AgilentEuroFemale, internal = TRUE, overwrite = TRUE)
    data(sysdata, envir=environment()) #load internal data
    if(reference != "male"){
      reference = AgilentEuroMale
    } else {
      reference = AgilentEuroFemale
    }
    
    names(reference)[6] = "GENOTYPE"
    names(reference)[7] = "DOUBLY_CUT"
    
    NonCurrent = which(data$genes$ProbeName %in% reference$PROBE_ID)
    current =  reference[-NonCurrent,-2]
    
    UncutCut = current[which(current$DOUBLY_CUT == 0),]
    UNCUT = which(paste(UncutCut$UNCUT_ALLELE,UncutCut$UNCUT_ALLELE, sep="") == as.character(UncutCut$GENOTYPE))
    probes = which(data$genes$ProbeName %in% UncutCut$PROBE_ID[UNCUT])
    
    if(array == "all") {
      dataSNP = data.frame(ID = 1:nrow(data$M[probes,]), data$M[probes,], stringsAsFactors = FALSE)
      dataSNP = melt(dataSNP, na.rm = TRUE, id = "ID" )
      names(dataSNP) = c("ID", "samples", "M")
    } else {
      dataSNP = data.frame(M = data$M[probes,array], stringsAsFactors = FALSE)
    }
  } else if(array != "all")
      dataSNP = data.frame(M = data$M[,array], stringsAsFactors = FALSE)
  
  type = match.arg(type)
  switch(type, 
         "log" = {
           if(array == "all"){
             ggplot(na.omit(dataSNP), aes(x=M, group = samples, colour = samples)) + 
               geom_density(adjust=1/5, fill = "transparent") +
               geom_vline(xintercept =mediana, lty = 2) +
               xlab("Log10 Ratio") + scale_y_sqrt() +
               ggtitle(mytitle)
           } else {
             ggplot(na.omit(dataSNP), aes(x=M)) + 
               geom_density(adjust=1/5, fill = brewer.pal(8,name = "Set2")[3]) +
               geom_vline(xintercept =mediana, lty = 2) +
               xlab("Log10 Ratio") + scale_y_sqrt() +
               ggtitle(mytitle)
           }
         },
         "theta" = {
           if(array == "all"){
             dataSNP$M = atan(dataSNP$M)
             ggplot(na.omit(dataSNP), aes(x=M, group = samples, colour = samples)) + 
               geom_density(adjust=1/5) +
               geom_vline(xintercept =mediana, lty = 2) +
               xlab("Theta") + scale_y_sqrt() +
               ggtitle(mytitle)
           } else {
             ggplot(atan(na.omit(dataSNP)), aes(x=M)) + 
               geom_density(adjust=1/5, fill = brewer.pal(8,name = "Set2")[2]) +
               geom_vline(xintercept =mediana, lty = 2) + 
               xlab("Theta") + scale_y_sqrt() +
               ggtitle(mytitle)
           }
           
         })
  
}