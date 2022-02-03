#tests which compare outputs of demofiles

test_that('Test demo.human1.R',{

    out=microPopModel(
        microbeNames=c('Bacteroides','NoButyStarchDeg','Acetogens'),
        times=seq(0,4,1/24),
        resourceSysInfo=resourceSysInfoHuman,
        microbeSysInfo=microbeSysInfoHuman,
        plotOptions=list(yLabel='concentration (g/l)',xLabel='time (d)',
            plotFig=FALSE,sumOverStrains=FALSE,saveFig=FALSE,
            figType='eps',figName='Human1')
    )

    finalrow=out$solution[nrow(out$solution),]

    expected.values=c(4.00,  3.44,  0.06,  0.03,  0.00,  0.00,  0.00,  3.20,  1.51,  2.40,  0.00, 10.01, 2.17, 10.34,  0.00,  0.00,    NA)

    expected.names=c('time','Bacteroides','NoButyStarchDeg','Acetogens','Protein','NSP','RS','Acetate','Propionate','Succinate','H2','CO2','other','H2O','Sugars','Formate','pH') 

    expect_equal(unname(round(finalrow,2)),expected.values)
    expect_equal(names(finalrow),expected.names)


})


test_that('Test demo.human2.R',{

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
        plotOptions=list(yLabel='concentration (g/l)',xLabel='time (d)',
            plotFig=FALSE,sumOverStrains=FALSE,saveFig=FALSE,
            figType='eps',figName='Human2'),
        pHLimit=TRUE
)


    finalrow=out$solution[nrow(out$solution),]

    
    expected.values=c( 4.00,  2.22,  1.24,  0.10,  0.00,  0.00 , 0.00,  4.86,  0.99 , 1.56 , 0.00, 10.01, 2.17, 10.01,  0.00,  0.00,  6.50)

    expected.names=c('time','Bacteroides','NoButyStarchDeg','Acetogens','Protein','NSP','RS','Acetate','Propionate','Succinate','H2','CO2','other','H2O','Sugars','Formate','pH') 

    expect_equal(unname(round(finalrow,2)),expected.values)
    expect_equal(names(finalrow),expected.names)


})

test_that('Test demo.human4.R',{


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
        numStrains=5,
        rateFuncs=myRateFuncs,
        plotOptions=list(yLabel='concentration (g/l)',xLabel='time (d)',
            plotFig=FALSE,sumOverStrains=FALSE,
            saveFig=FALSE,figType='eps',figName='Human4'),
        pHLimit=TRUE,
        strainOptions=list(
            randomParams=c('halfSat','yield','maxGrowthRate','pHtrait'),
            seed=1,
            distribution='uniform',
            percentTraitRange=5,
            maxPHshift=0.05,
            applyTradeOffs=TRUE,
            tradeOffParams=c('halfSat','maxGrowthRate'),
            paramsSpecified=TRUE,
            paramDataName=strainParams)
    ) 

    finalrow=out$solution[nrow(out$solution),]

    expected.names=c('time',paste0('Bacteroides.',seq(1,5)),paste0('NoButyStarchDeg.',seq(1,5)),paste0('Acetogens.',seq(1,5)),'Protein','NSP','RS','Acetate','Propionate','Succinate','H2','CO2','other','H2O','Sugars','Formate','pH') 

    expected.values=c(4.00,  0.11 , 0.09,  0.09 , 0.08,  1.95,  0.13 , 0.62,  0.06,  0.30 , 0.07,  0.03,  0.02 , 0.02,  0.02 , 0.02,  0.00,  0.00 , 0.00 , 4.80,  1.00 , 1.58,  0.00, 10.00,  2.15 ,10.02 , 0.00,  0.00,  6.50)

    expect_equal(unname(round(finalrow,2)),expected.values)
    expect_equal(names(finalrow),expected.names)

})
