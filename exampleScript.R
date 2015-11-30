########################
#' @_Example_Script
#' 
#' @_DO_NOT_RUN
########################

install.packages("ggdendro")

setwd("D:/Sources/Skryabin")
RG = readAgilent(makeTargets()[1:4,], format = "F")

setwd("D:/Sources/Vasilyev/may2015F")
E = readAgilent(makeTargets()[1:4,], format = "F", dyes = 1)

#' @_PLOT_RGList

gplotMA(RG)
gplotDendro(RG)
gplotKMeans(RG)
gplotDensities(RG)

#' @_PLOT_RGList_with_SNP

gplotSNP(RG)
gplotSNP(RG, type = "theta")
gplotSNP(RG, array = "all", type = "theta")
gplotSNP(RG, type = "theta", SNPprobe = 1)

#' @_PLOT_RGList_with_SNP
RG$genes$ProbeName[7:20] = c("A_20_P00100012", "A_20_P00201918", "A_20_P00100018", "A_20_P00201924", "A_20_P00201926", "A_20_P00100020",
                             "A_20_P00201929", "A_20_P00201931", "A_20_P00201932", "A_20_P00100026", "A_20_P00201933", "A_20_P00201936",
                             "A_20_P00100030", "A_20_P00201938")
gplotSNP(RG, reference = "male", SNPprobe = 0, type = "theta")

#' @_PLOT_RGList

gplotDendro(E)
gplotKMeans(E)
gplotDensities(E)

#' @_outliners

findOutliners(RG)
findOutliners(RG, flags = FALSE)
gplotOutliners(RG, 3)
gplotOutliners(RG, 1, FALSE)

RG2 = rmOutliners(RG) 

#

findOutliners(E)
findOutliners(E, flags = FALSE)
gplotOutliners(E, 3)

E2 = rmOutliners(E, method = "between")
gplotOutliners(E2, flags = FALSE)

#

gplotArray(RG2, 1)
gplotArray(E2, 3)

#' @_spatial_correction

gplotRC(RG)
gplotRC(E)

RGf = fixArrays(RG)
Ef = fixArrays(E)

gplotRC(RG,1)
gplotRC(Ef)

################################
CNVnorm = MA.RG(RGf)
RGf2 = noiseReductionCGH(RGf)
boxplot(log(RGf2$R), col="red")
boxplot(log(RGf2$G), add=T, col="green")

CNVnorm = normalizeWithinArrays(RGf2, method = "loess", bc.method = "none")
CNVnorm2 = normalizeBetweenArrays(CNVnorm)

boxplot(MA.RG(RG)$M)
boxplot(MA.RG(RGf2)$M, add=T, col=2)
boxplot(CNVnorm$M, add=T, col=3)

gplotMA(RG,1)
gplotMA(RGf,1)
gplotMA(RGf2,1)

##
RGc = RGf[RG$genes$ControlType == 1,]
gplotMA(RGc,1)
##

gplotMA(CNVnorm,1)
gplotMA(CNVnorm2,1)
ma.plot(na.omit(cbind(CNVnorm$A[,smp], CNVnorm$M[,smp]))[,1],
        na.omit(cbind(CNVnorm$A[,smp], CNVnorm$M[,smp]))[,2],plot.method = "smoothScatter")


gplotDensities(MA.RG(RG))
gplotDensities(MA.RG(RGf))
gplotDensities(CNVnorm)
gplotDensities(CNVnorm2)

gplotRC(RG.MA(CNVnorm),1)
gplotRC(RG.MA(CNVnorm2),1)

Enorm = nec(Ef, status = Ef$genes$ControlType, negctrl="-1",  regular="0",  offset=1, robust=TRUE)
Enorm2 = normalizeBetweenArrays(Enorm, method="quantile")

gplotDensities(E)
gplotDensities(Ef)
gplotDensities(Enorm)
gplotDensities(Enorm2)

gplotRC(Enorm,4)
gplotRC(Enorm2,4)


####### array CGH

DLRSpread(CNVnorm)
cna = runCBS(CNVnorm, cluster = 2)
cnas = thrCBS(cna, n=0.5)
statCNA(cnas)
csd = commonSeg(cnas, cna, biomart = TRUE)
data.frame(apply(csd[1:10,], 2, strtrim, 8))
# write.table(CSD, file = "CSD.txt", quote = F, sep = "\t", na = "", row.names = F)

######## plots CGH

setwd("D:/Sources/Skryabin/genomes")
ggenomeplot(cna, smp = 1, span = 0.3)


######## CGh call (work with tumor cells)

CGHc = convertCGHcall(CNVnorm)
head(CGHc);tail(CGHc)
library(CGHcall)
# average log2 transformation in duplicate probe
# find colnames that not need averaging process
# col.not.need <- colnames(CGHc)[!grepl('kap_7',colnames(CGHc))]
#convert data.frame to data.table
#library(data.table)
#tmp.callRaw = as.data.table(CGHc)
#averaging
# ave.callRaw = tmp.callRaw[,list('kap_7'= median(kap_7)),col.not.need]
# ave.callRaw = tmp.callRaw[,lapply(.SD,median),col.not.need] 

callRaw.maked <- make_cghRaw(CGHc)
call.cghdata <- preprocess(callRaw.maked, maxmiss=30, nchrom=22)
call.norm.cghdata <- normalize(call.cghdata, method="median", smoothOutliers=TRUE)
seg.calldata <- segmentData(call.norm.cghdata, method="DNAcopy",
                            undo.splits="sdundo",undo.SD=3,clen=10, relSDlong=5)

postseg.calldata <- postsegnormalize(seg.calldata)

cellularity = c(1, 1) #threshold
result.calling <- CGHcall(postseg.calldata,nclass=3,cellularity=cellularity)
# convert to a call object
result.calling <- ExpandCGHcall(result.calling,postseg.calldata)

#images
plot(result.calling[,1])
plot(result.calling[,2])
plot(result.calling[,3])
plot(result.calling[,4])
summaryPlot(result.calling)
frequencyPlotCalls(result.calling)



