rm(list=ls())
graphics.off()

detach(package:microPop)
library(microPop)


out=microPopModel(
    microbeNames=c('Bacteroides','ButyrateProducers1'),
    times=seq(0,10,0.1),
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeSysInfoHuman,
    numStrains=c(Bacteroides=2,ButyrateProducers1=4),
    strainOptions=list(randomParams=c('maxGrowthRate','halfSat'),percentTraitRange=10),
    plotOptions=list(sumOverStrains=TRUE,plotFig=TRUE),
    pHLimit=TRUE,
    pHVal=5.7
)


