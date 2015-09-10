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

gplotRC(RG,4)
gplotRC(E)

RGf = fixArrays(RG)
RGf = fixArrays(RGf)
Ef = fixArrays(E)

gplotRC(RGf,4)
gplotRC(Ef)

################################

CNVnorm = normalizeWithinArrays(RGf, method = "loess", iterations = 20, bc.method = "none")
CNVnorm2 = normalizeBetweenArrays(CNVnorm)

gplotMA(RG,1)
gplotMA(RGf,2)
gplotMA(CNVnorm,2)
gplotMA(CNVnorm2,2)

gplotDensities(MA.RG(RG))
gplotDensities(MA.RG(RGf))
gplotDensities(CNVnorm)
gplotDensities(CNVnorm2)

gplotRC(CNVnorm,4)
gplotRC(CNVnorm2,4)

Enorm = nec(Ef, status = Ef$genes$ControlType, negctrl="-1",  regular="0",  offset=1, robust=TRUE)
Enorm2 = normalizeBetweenArrays(Enorm, method="quantile")

gplotDensities(E)
gplotDensities(Ef)
gplotDensities(Enorm)
gplotDensities(Enorm2)

gplotRC(Enorm,4)
gplotRC(Enorm2,4)


