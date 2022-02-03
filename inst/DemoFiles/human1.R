#human colon 1
#Simplest example - look at 3 MFGs: 'Bacteroides','NoButyStarchDeg','Acetogens'
#No dependence on pH


library(microPop)


out=microPopModel(
    microbeNames=c('Bacteroides','NoButyStarchDeg','Acetogens'),
    times=seq(0,4,1/24),
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeSysInfoHuman,
    plotOptions=list(yLabel='concentration (g/l)',xLabel='time (d)',
                     plotFig=TRUE,sumOverStrains=FALSE,saveFig=FALSE,
                     figType='eps',figName='Human1',cex.plot=1,cex.legend=0.7)
)


print(out$solution[1,])

