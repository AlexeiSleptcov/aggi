#' Calculate common CNA statistic 
#'
#' Calculate common CNA statistic 
#' 
#' @param data CBS data, found CNA from \code{thrCBS}
#' @param data2 DNAcopy data
#' @param biomart logical, if TRUE \code{\link{biomaRt}} will be used for description
#'
#' @export

commonSeg <- function(data, data2, biomart=FALSE){
  if(!inherits(data, "CBS"))
    stop("First, data must be CBS", call.=0)
  if(!inherits(data2, "DNAcopy"))
    stop("First, data must be CBS", call.=0)
  if(exists("data3"))
    if(!inherits(data3, c("MAList", "RGList")))
      stop("First, data must be CBS", call.=0)
  
  
  dl = do.call(rbind, data)
  # sorting data by maploc
  dl = dl[order(dl$loc.end),]
  dl = dl[order(dl$loc.start),]
  dl = dl[order(dl$chrom),]
  rownames(dl) = 1:nrow(dl)
  
  td = data.frame(rbind(rep(NA,length(data))))
  colnames(td) = names(data)
  dd = cbind(dl, td, "sd"=0, td)[0,]
  
  for(chr in c(1:22,"X","Y")){
    if(!sum(dl$chrom == chr)) next;
    da = dl[dl$chrom == chr,]
    d2 = data2$data[data2$data$chrom == chr,]
    
    # iter by row
    i = 1
    repeat{
      if(i > nrow(da)) 
        break;
      da.s = da$loc.start %in% da$loc.start[i]
      da.e = da$loc.end %in% da$loc.end[i]
      who  = which(da.s & da.e)
      nam  = da[who,"ID"]
      td2  = td
      td2[colnames(td2) %in% nam] = round(da[who,"seg.mean"], 2)
      # where maploc
      w.ml = suppressWarnings(which(d2$maploc == da[i,"loc.start"]):which(d2$maploc == da[i,"loc.end"]))
      ts   = t(data.frame("X1" = round(apply(d2[w.ml, 3:ncol(d2)], 2, median, na.rm=1), 2)))
      # save result
      dd2 = cbind(da[i,], td2, "sd"=0, ts)
      dd2$sd = round(sd(da[who,"seg.mean"]), 2)
      dd  = rbind(dd, dd2)
      
      # mechanic
      if(length(who) > 1) {
        i = i+length(who)
      } else {
        i = i+1
      }
    }# end repeat
  }#end for
  
  dd_1 = dd[,2:4]
  dd_2 = dd[,5:ncol(dd)]
  colnames(dd_1) = c("chr", "start", "end")
  colnames(dd_2)[1] = "N"
  # get size
  dd_1 = cbind(dd_1, "size" = dd_1$end-dd_1$start+1)
  
  # describe maploc
  if(biomart){
    require("biomaRt")
    attr = c('band','entrezgene', 'hgnc_id', 'hgnc_symbol')
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    filters = listFilters(ensembl)[1:3,1]
    dd_1_2 = cbind(dd_1,matrix(NA,nrow=nrow(dd_1),ncol=length(attr),dimnames=list(1:nrow(dd_1),attr)))
    
    # create progessBar
    pb <- txtProgressBar(min = 0, max = nrow(dd_1_2), style = 3, width=50) # progress bar
    
    for(i in 1:nrow(dd_1_2)){
      # request from dd_1_2
      request.value = list(dd_1_2[i,'chr'],dd_1_2[i,'start'],dd_1_2[i,'end'])
      ans.BM = getBM(attributes=attr, filters=filters, values=request.value, mart=ensembl, bmHeader=T)
      # ckeck result
      if(nrow(ans.BM)){
        
        # check wide band & knock into one
        b = unique(ans.BM$Band)
        if(length(b) == 1){
          dd_1_2$band[i] = b
        } else {
          dd_1_2$band[i] = paste(b[1], b[length(b)], sep="-")
        }
        # check multiple attr & knock into one
        for(k in 2:length(attr))
          dd_1_2[[attr[k]]][i] = paste(unique(na.omit(ans.BM[[k]])), collapse=", ")     
      } 
      setTxtProgressBar(pb, i) # set progressBar
      # delete "hgnc_id:"
      dd_1_2$hgnc_id = gsub("HGNC:([1:9]*)", "\\1", dd_1_2$hgnc_id) 
      res = cbind(dd_1_2, dd_2)
    } #end for
    close(pb) # end progessBar 
    ###   need data3 (e.g. BList) ###
    #   } else if(exists("data3")){
    #     data3  = data3[data3$genes$ControlType == 0,]
    #     starts = suppressWarnings(maploc(data3))
    #     chr    = chroma(data3)
    #     genes  = c()
    #     for(i in 1:nrow(dd)){
    #       w.st  = suppressWarnings(which(chr %in% dd[i,"chrom"] & starts %in% dd[i,"loc.start"]:dd[i,"loc.end"]))
    #       gene1 = paste(unique(data3$genes$GeneName[w.st]), collapse = ", ")
    #       genes = rbind(genes, gene1)
    #       res = cbind(dd_1, genes, dd_2)
    #     }
  }# end desc
  res
}
