## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(microPop)

## ----echo=FALSE---------------------------------------------------------------
#install.packages('~/MicroPop/MicroPop1/microPop_1.6.tar.gz',repos=NULL,source=TRUE)
library(knitr)

## -----------------------------------------------------------------------------
library(microPop)

## -----------------------------------------------------------------------------
#check package version 
packageVersion("microPop")

## -----------------------------------------------------------------------------
M1=createDF(file='MFG1.csv')

## ----echo=FALSE---------------------------------------------------------------
kable(M1)

## ---- echo=FALSE, fig.height=4,fig.width=4------------------------------------
s=seq(0,0.05,0.0001)
Gmax=5
K=0.001
G=Gmax*s/(s+K)
  plot(range(s),c(0,5),type='n',xlab='substrate (g/l)',ylab='Specific growth rate (/d)')
  lines(s,G)
  abline(v=K,lty=2,col=2)
  abline(h=Gmax,col=3,lty=2)
  legend('right',col=c(2,3),lty=2,legend=c('halfSat','maxGrowthRate'),bty='n')

## -----------------------------------------------------------------------------
res.SI=createDF(file='resourceSI.csv')

## ----echo=FALSE---------------------------------------------------------------
kable(res.SI)

## -----------------------------------------------------------------------------
mic.SI=createDF(file="microbeSI.csv")

## ---- echo=FALSE--------------------------------------------------------------
kable(mic.SI)

## ---- results='hide'----------------------------------------------------------
out=microPopModel(microbeNames = 'M1',
              times=seq(0,1,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=mic.SI,
              plotOptions=list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8------------------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## -----------------------------------------------------------------------------
names(out)

## -----------------------------------------------------------------------------
head(out$solution)

## ---- fig.width=5, fig.height=4-----------------------------------------------
plot(range(out$solution[,'time']),c(0,10),type='n',xlab='Time (d)',ylab='concentration (g/l)')
lines(out$solution[,'time'],out$solution[,'M1'],col='black')
lines(out$solution[,'time'],out$solution[,'S1'],col='red')
lines(out$solution[,'time'],out$solution[,'P1'],col='blue')
legend('topright',c('M1','S1','P1'),col=c('black','red','blue'),lty=1)

## -----------------------------------------------------------------------------
names(out$parms)

## -----------------------------------------------------------------------------
out$parms$microbeNames
out$parms$resourceNames

## -----------------------------------------------------------------------------
M2=createDF(file="MFG2.csv")

## ----echo=FALSE---------------------------------------------------------------
kable(M2)

## -----------------------------------------------------------------------------
micNew.SI = cbind(mic.SI,M2=mic.SI[,'M1']) #cbind is a function which binds columns

## ----echo=FALSE---------------------------------------------------------------
kable(micNew.SI)

## ---- results='hide'----------------------------------------------------------
out=microPopModel(microbeNames = c('M1','M2'),
              times=seq(0,3,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=micNew.SI,
              plotOptions = list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8, echo=FALSE------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION 
#look at the final concentrations:

out$solution[nrow(out$solution),c('M1','M2')]

#M1 and M2 have the same yield but M2 has a higher max growth rate and a higher halfSat. 
#As time goes on the substrate concentration is low so the max growth rate is not important.
#We change half-sat of M2 to 1e-4 (i.e. less than M1) and look what happens:
newM2=M2 #make a new copy of M2 to play with
newM2['halfSat','S1']=1e-4
micNew.SI = cbind(micNew.SI,newM2=micNew.SI[,'M2']) #add newM2 to the microbe SI data frame (assume same conditions as M2)

out=microPopModel(microbeNames = c('M1','newM2'),
                  times=seq(0,3,0.001),
                  resourceSysInfo=res.SI,
                  microbeSysInfo=micNew.SI,
                  plotOptions = list(plotFig=FALSE)
)

out$solution[nrow(out$solution),c('M1','newM2')]
#Now M2 > M1 near the end of the simulation

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION 
out=microPopModel(microbeNames = c('M1','M2'),
                  times=seq(0,10,0.001),
                  resourceSysInfo=res.SI,
                  microbeSysInfo=micNew.SI,
                  plotOptions = list(plotFig=FALSE)
)

out$solution[nrow(out$solution),c('M1','M2')]
#M2 dies out as it can't compete with M1 at low substrate concentrations.

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION
#start by making M1 and M2 the same:
newM2[c('halfSat','yield','maxGrowthRate'),'S1'] = M1[c('halfSat','yield','maxGrowthRate'),'S1'] 

out=microPopModel(microbeNames = c('M1','newM2'),
                  times=seq(0,6,0.001),
                  resourceSysInfo=res.SI,
                  microbeSysInfo=micNew.SI,
                  plotOptions=list(plotFig=FALSE)
)

#Look at steady state values
out$solution[nrow(out$solution),c('M1','newM2','S1')]


#Now investigate the effect of yield values:

#make M2 have a lower yield than M1
newM2['yield','S1'] = 0.5*as.numeric(M1['yield','S1'])

out=microPopModel(microbeNames = c('M1','newM2'),
                  times=seq(0,6,0.001),
                  resourceSysInfo=res.SI,
                  microbeSysInfo=micNew.SI,
                  plotOptions = list(plotFig=FALSE)
)
out$solution[nrow(out$solution),c('M1','newM2','S1')]

## ---- results='hide'----------------------------------------------------------
out=microPopModel(microbeNames = c('M1','M2'),
              times=seq(0,3,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=micNew.SI,
              numStrains=3,
              strainOptions=list(randomParams=c('halfSat','yield','maxGrowthRate'),
                                 seed=1,distribution='uniform',
                                 percentTraitRange=20,
                                 applyTraitTradeOffs=TRUE,
                                 tradeOffParams=c('halfSat','maxGrowthRate')),
              plotOptions = list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8------------------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7,sumOverStrains = FALSE)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## ---- fig.width=5, fig.height=4-----------------------------------------------
barplot(out$solution[nrow(out$solution),out$parms$allStrainNames],col=c(rep('red',3),rep('cyan',3)))

## -----------------------------------------------------------------------------
out$parms$allStrainNames

## ---- results='hide'----------------------------------------------------------
out=microPopModel(microbeNames = c('M1','M2'),
              times=seq(0,3,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=micNew.SI,
              numStrains=c(M1=2,M2=3),
              strainOptions=list(randomParams=c('halfSat','yield','maxGrowthRate'),
                                 seed=1,distribution='uniform',
                                 percentTraitRange=20,
                                 applyTraitTradeOffs=TRUE,
                                 tradeOffParams=c('halfSat','maxGrowthRate')),
              plotOptions = list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8,echo=FALSE-------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7,sumOverStrains = FALSE)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION
out3=out #save the 3 day solution

out10=microPopModel(microbeNames = c('M1','M2'),
              times=seq(0,10,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=micNew.SI,
              numStrains=3,
              strainOptions=list(randomParams=c('halfSat','yield','maxGrowthRate'),
                                 seed=1,distribution='uniform',
                                 percentTraitRange=20,
                                 applyTraitTradeOffs=TRUE,
                                 tradeOffParams=c('halfSat','maxGrowthRate')),
              plotOptions = list(plotFig=FALSE)
)


## ----echo=TRUE, eval=TRUE, fig.width=5, fig.height=4--------------------------
barplot(out10$solution[nrow(out10$solution),out10$parms$allStrainNames],col=c(rep('red',3),rep('cyan',3)))

## ----echo=TRUE, eval=TRUE-----------------------------------------------------

out20=microPopModel(microbeNames = c('M1','M2'),
              times=seq(0,20,0.001),
              resourceSysInfo=res.SI,
              microbeSysInfo=micNew.SI,
              numStrains=3,
              strainOptions=list(randomParams=c('halfSat','yield','maxGrowthRate'),
                                 seed=1,distribution='uniform',
                                 percentTraitRange=20,
                                 applyTraitTradeOffs=TRUE,
                                 tradeOffParams=c('halfSat','maxGrowthRate')),
              plotOptions = list(plotFig=FALSE)
)


## ----echo=TRUE, eval=TRUE, fig.width=5, fig.height=4--------------------------
barplot(out20$solution[nrow(out20$solution),out20$parms$allStrainNames],col=c(rep('red',3),rep('cyan',3)))

## ----echo=TRUE, eval=TRUE, fig.width=5, fig.height=4--------------------------
#Plot all together
barplot(rbind(out3$solution[nrow(out3$solution),out3$parms$allStrainNames],
              out10$solution[nrow(out10$solution),out10$parms$allStrainNames],
              out20$solution[nrow(out20$solution),out20$parms$allStrainNames]),
        beside=TRUE,legend.text = c('3 days','10 days','20 days'))


## -----------------------------------------------------------------------------
  names(rateFuncsDefault)

## -----------------------------------------------------------------------------
  rateFuncsDefault$growthLimFunc

## ---- eval=FALSE--------------------------------------------------------------
#    myRateFuncs = rateFuncsDefault #make your own copy of the rateFuncs list
#    myRateFuncs$growthLimFunc = myGrowthFunction #assign your own function

## ---- eval=FALSE--------------------------------------------------------------
#  out=microPopModel(microbeNames = 'M1',
#                times=seq(0,10,0.001),
#                resourceSysInfo=res.SI,
#                microbeSysInfo=mic.SI,
#                rateFuncs=myRateFuncs
#  )

## -----------------------------------------------------------------------------
M3=createDF(file='MFG3.csv')

## ----echo=FALSE---------------------------------------------------------------
kable(M3)

## ----fig.width=5,fig.height=3-------------------------------------------------
plot(as.numeric(M3['pHcorners',2:5]),c(0,1,1,0),type='l',xlab='pH',ylab='pH limitation',col='blue',lwd=2)

## -----------------------------------------------------------------------------
res.SI2=createDF(file='resourceSI2.csv')

## ----echo=FALSE---------------------------------------------------------------
kable(res.SI2) 

## -----------------------------------------------------------------------------
mic.SI2=createDF(file='microbeSI2.csv')

## ----echo=FALSE---------------------------------------------------------------
kable(mic.SI2) 

## ---- results='hide', fig.show='hide'-----------------------------------------
out1=microPopModel(microbeNames = 'M3',
              times=seq(0,5,0.001),
              resourceSysInfo=res.SI2,
              microbeSysInfo=mic.SI2,
              pHLimit=TRUE,
              pHVal=4.1,
              plotOptions=list(plotFig=FALSE)
              )
    
out2=microPopModel(microbeNames = 'M3',
              times=seq(0,5,0.001),
              resourceSysInfo=res.SI2,
              microbeSysInfo=mic.SI2,
              pHLimit=TRUE,
              pHVal=6,
              plotOptions=list(plotFig=FALSE)
              )

## ---- fig.width=5,fig.height=4------------------------------------------------
plot(range(out1$solution[,'time']),c(0,max(c(out1$solution[,'M3'],out2$solution[,'M3']))),type='n',xlab='Time (d)',ylab='Concentration (g/l)', main='M3',lwd=2)
lines(out1$solution[,'time'],out1$solution[,'M3'],col='black',lwd=2)
lines(out2$solution[,'time'],out2$solution[,'M3'],col='red',lwd=2)
legend('right',col=1:2,lty=1,legend=c(4.1,6),title='pH',lwd=2)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION
myRateFuncs=rateFuncsDefault

myRateFuncs$pHFunc=function(time,parms){
  if (time<=2){
    pH=6
  }else{
    pH=4.1
  }
  return(pH)
}

out=microPopModel(microbeNames = 'M3',
              times=seq(0,5,0.001),
              resourceSysInfo=res.SI2,
              microbeSysInfo=mic.SI2,
              pHLimit=TRUE,
              rateFuncs=myRateFuncs
)

## -----------------------------------------------------------------------------
data(package='microPop')

## ----results='hide'-----------------------------------------------------------
system.file('DemoFiles/',package='microPop')

## ---- eval=FALSE--------------------------------------------------------------
#  runMicroPopExample('human1')

## ---- results='hide', fig.show='hide'-----------------------------------------
out=microPopModel(
    microbeNames=c('Bacteroides','NoButyStarchDeg','Acetogens'),
    times=seq(0,4,1/24),
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeSysInfoHuman,
    plotOptions = list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8,echo=FALSE-------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## ---- results='hide',fig.show='hide'------------------------------------------

simulation.times=seq(0,4,1/24)

myRateFuncs=rateFuncsDefault

myRateFuncs$pHFunc=function(time,parms){
    if (time<=max(simulation.times)/2){pH=5.5}else{pH=6.5}
    return(pH)
}

out=microPopModel(
    microbeNames=c('Bacteroides','NoButyStarchDeg','Acetogens'),
    times=simulation.times,
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeSysInfoHuman,
    rateFuncs=myRateFuncs,
    pHLimit=TRUE,
    plotOptions = list(plotFig=FALSE)
)

## ---- fig.height=4,fig.width=8,echo=FALSE-------------------------------------
par(mfrow=c(1,2))
plotMicrobes(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)
plotResources(out,cex.title=0.7,cex.legend=0.5,cex.ax=0.7)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
#SOLUTION

myBacteroides=Bacteroides

print(myBacteroides['pHcorners',])
print(myBacteroides['pHcorners',2:5])
myBacteroides['pHcorners',2:5]=c(5,6,7,8)

microbeHuman.SI = cbind(microbeSysInfoHuman,myBacteroides=microbeSysInfoHuman[,'Bacteroides']) #add Sys Info

out=microPopModel(
    microbeNames=c('myBacteroides','NoButyStarchDeg','Acetogens'),
    times=simulation.times,
    resourceSysInfo=resourceSysInfoHuman,
    microbeSysInfo=microbeHuman.SI,
    rateFuncs=myRateFuncs,
    pHLimit=TRUE
)


## ---- echo=FALSE--------------------------------------------------------------
kable(M3)

## ---- echo=FALSE--------------------------------------------------------------
kable(Xh2)

