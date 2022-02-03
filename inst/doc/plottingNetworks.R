## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----results='hide'-----------------------------------------------------------
library(microPop)
library(visNetwork)

## -----------------------------------------------------------------------------
v1=networkDFfromMPinput('Bacteroides') #this produces the network nodes and edges dataframes (contained in the list v1)
v2=getVNPlotObject(v1$nodes,v1$edges) #this takes the nodes and edges and makes the visNetwork object for plotting
visOptions(v2) #this produces a plot of the network. You can also use print(v2)

## -----------------------------------------------------------------------------
data(package='microPop')

## -----------------------------------------------------------------------------
microbes=c('Acetogens','ButyrateProducers1','Bacteroides','Methanogens')
DFs=networkDFfromMPinput(microbes)
nodes=DFs$nodes
edges=DFs$edges
vv=getVNPlotObject(nodes,edges,
                    addLegend=TRUE,
                    figWidth=600,
                    figHeight=600)

visOptions(vv)

## -----------------------------------------------------------------------------
vv=getVNPlotObject(nodes,edges,
                   addLegend=TRUE,
                   microbeCol='pink',
                   resourceCol='cyan',
                   productionCol='red',
                   uptakeCol='purple',
                   figWidth=600,
                   figHeight=600)

visOptions(vv)

## -----------------------------------------------------------------------------
allMicrobes=c('Acetogens','ButyrateProducers1','ButyrateProducers2','ButyrateProducers3',
                           'Bacteroides','Methanogens','LactateProducers',
              'NoButyFibreDeg','NoButyStarchDeg','PropionateProducers')

DFs=networkDFfromMPinput(allMicrobes)
nodes=DFs$nodes
edges=DFs$edges

#get the visNetwork plot object:
vv=getVNPlotObject(nodes,edges,
                   addLegend=FALSE,
                   mainTitle='10 groups',
                   figWidth=800,
                   figHeight=800)

visOptions(vv)

## ---- results='hide'----------------------------------------------------------

  microbes=c('Acetogens','Bacteroides','ButyrateProducers1')

  out=microPopModel(microbeNames=microbes,
                      times=seq(0,2,0.01),
                      resourceSysInfo=resourceSysInfoHuman,
                      microbeSysInfo=microbeSysInfoHuman,
                      plotOptions = list(plotFig=FALSE),
                      networkAnalysis=TRUE)

## -----------------------------------------------------------------------------
DFs=networkDFfromMPoutput(chosen.time=1,MPoutput=out)

## -----------------------------------------------------------------------------
vv=getVNPlotObject(DFs$nodes,DFs$edges,scaleNodes=TRUE,scaleEdges = TRUE,
                   mainTitle='Network after 1 d', 
                   figWidth = 500,figHeight = 500)
visOptions(vv)

## ---- results='hide'----------------------------------------------------------

  microbes=c('Acetogens','Bacteroides','ButyrateProducers1')

  out=microPopModel(microbeNames=microbes,
                      times=seq(0,2,0.01),
                      numStrains=c(Acetogens=1,Bacteroides=2,ButyrateProducers1=3),
                      strainOptions=list(randomParams=c('maxGrowthRate','halfSat'),percentTraitRange=20),
                      resourceSysInfo=resourceSysInfoHuman,
                      microbeSysInfo=microbeSysInfoHuman,
                      plotOptions = list(plotFig=FALSE),
                      networkAnalysis=TRUE)

## -----------------------------------------------------------------------------
DFs=networkDFfromMPoutput(chosen.time=1,MPoutput=out,sumOverStrains = FALSE)

vv=getVNPlotObject(DFs$nodes,DFs$edges,scaleNodes=TRUE,scaleEdges = TRUE,
                   figWidth = 700,figHeight = 700)
visOptions(vv)

