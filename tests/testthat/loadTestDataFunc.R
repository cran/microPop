loadTestDataFunc=function(testCase){

    library(microPop)
    
    #only run microPopModel() when a different system is required
    #Note test.flag is undefined until microPopModel() has been called

    #testCase=0
    #any type of microPop run

    #testCase=1
    #essential resources and only one strain

    #testCase=2
    #multiple strains

    #testCase=3
    #simple case for comparison with analytical solution

    #testCase=4
    #groups with differing numbers of strains
    
 
    
    if (testCase==1 | testCase==0){ #1 strain, essential resources-------------------------
        
        MFG=matrix(NA,ncol=4,nrow=6,dimnames=list(c('Rtype','halfSat','yield',
                                        'maxGrowthRate','stoichiom','keyResource'),
                                        c('H2','CO2','CH4','H2O')))
        MFG['Rtype',]=c('Se','Se','P','P')
        MFG['halfSat',c('H2','CO2')]=1e-6
        MFG['yield','H2']=0.2
        MFG['maxGrowthRate','H2']=2
        MFG['keyResource',1]='H2'
        MFG['stoichiom',]=c(4,1,1,2)
        Archea <<- data.frame(MFG,stringsAsFactors=FALSE)
        
        Rmat=matrix(NA,ncol=4,nrow=4,dimnames=list(c('startValue','inflowRate',
                                         'washOut','molarMass'),c('H2','CO2','CH4','H2O')))
        Rmat['startValue',]=c(1,1,0,0)
        Rmat['inflowRate',]=c(1,5,0,0)
        Rmat['washOut',]=c(0.1,0.1,0.1,0.1)
        Rmat['molarMass',]=c(2,44,16,18)
        resInfo1 <<- data.frame(Rmat,stringsAsFactors=FALSE)
        
        Mmat=matrix(NA,ncol=1,nrow=3,dimnames=list(c('startValue','inflowRate',
                                         'washOut'),c('Archea')))
        Mmat['startValue',]=1
        Mmat['inflowRate',]=0
        Mmat['washOut',]=0.1
        micInfo1 <<- data.frame(Mmat,stringsAsFactors=FALSE)
        
        resInfo=resInfo1
        micInfo=micInfo1
        microbeNames='Archea'
        numStrains=1
        strainOptions=list(randomParams=NULL)
        pHLim=FALSE
        pHval=NULL
        times=seq(0,50,0.1)
    }

    
    if (testCase==2){ #multiple strains-------------------------------------------------
        resInfo=resourceSysInfoHuman
        micInfo=microbeSysInfoHuman
        microbeNames=c('Bacteroides','Acetogens')
        numStrains=3
        strainOptions<<-list(
            randomParams=c('halfSat','yield','maxGrowthRate','pHtrait'),
            seed=1,
            distribution='uniform',
            percentTraitRange=5,
            maxPHshift=0.05,
            applyTradeOffs=TRUE,
            tradeOffParams=c('halfSat','maxGrowthRate'))
        
        pHLim=TRUE
        pHval=6.5
        times=seq(0,50,0.1)
    }


    if (testCase==3){#simple case for comparison with analytical solution------------------
        Gmax=10
        K=0.01
        V=2
        Y=0.3
        SVin=10
        
        init=c(S=1,B=1,P=0)

        MFG=matrix(NA,ncol=2,nrow=5,dimnames=list(c('Rtype','halfSat','yield',
                                            'maxGrowthRate','stoichiom'),c('S','P')))
        MFG['Rtype',]=c('S','P')
        MFG[c('halfSat','yield','maxGrowthRate'),'S']=c(K,Y,Gmax)
        MFG['stoichiom',]=c(1,1)
        B<<-data.frame(MFG,stringsAsFactors=FALSE)
  
    #make resourceSysInfo data frame
        Rmat=matrix(NA,ncol=2,nrow=4,dimnames=list(c('startValue','inflowRate',
                                             'washOut','molarMass'),c('S','P')))
        Rmat['startValue',]=init[c('S','P')]
        Rmat['inflowRate',]=c(SVin,0)
        Rmat['washOut',]=V
        Rmat['molarMass',]=c(1,1)
        
        Mmat=matrix(NA,ncol=1,nrow=3,dimnames=list(c('startValue','inflowRate','washOut'),c('B')))
        Mmat['startValue',]=init['B']
        Mmat['inflowRate',]=0
        Mmat['washOut',]=V
        
        microbeNames='B'
        times=seq(0,10,0.001)
        resInfo<<-data.frame(Rmat,stringsAsFactors=FALSE)
        micInfo<<-data.frame(Mmat,stringsAsFactors=FALSE)
        numStrains=1
        strainOptions=list(randomParams=NULL)
        pHLim=FALSE
        pHval=NULL
    }

     if (testCase==4){ #same as case 2 but with different number of strains in each group
         resInfo=resourceSysInfoHuman
         micInfo=microbeSysInfoHuman
         microbeNames=c('Bacteroides','Acetogens')
         numStrains=c(Bacteroides=3,Acetogens=2)
         strainOptions<<-list(
             randomParams=c('halfSat','yield','maxGrowthRate','pHtrait'),
             seed=1,
             distribution='uniform',
             percentTraitRange=5,
             maxPHshift=0.05,
             applyTradeOffs=TRUE,
             tradeOffParams=c('halfSat','maxGrowthRate'))
         
         pHLim=TRUE
         pHval=6.5
         times=seq(0,50,0.1)
     }


    #---------------------------------------------------------------------------------

    
    if (exists('test.flag')){ #i.e. function has been run before

        if (testCase==0 & test.flag==2){
            resInfo=resourceSysInfoHuman
            micInfo=microbeSysInfoHuman
            microbeNames=c('Bacteroides','Acetogens')
            numStrains=3
            strainOptions<<-list(
                randomParams=c('halfSat','yield','maxGrowthRate','pHtrait'),
                seed=1,
                distribution='uniform',
                percentTraitRange=5,
                maxPHshift=0.05,
                applyTradeOffs=TRUE,
                tradeOffParams=c('halfSat','maxGrowthRate'))
            
            pHLim=TRUE
            pHval=6.5
        }


        #==========================run microPop=========================
        
        #don't need to rerun any thing if testCase is 0 as system type doesn't matter
        if (testCase != 0){
            
            if (test.flag != testCase){ #but not for right case
                
                out <<- microPopModel(
                    microbeNames=microbeNames,
                    times=times,
                    resourceSysInfo=resInfo,
                    microbeSysInfo=micInfo,
                    plotOptions=list(plotFig=FALSE),
                    numStrains=numStrains,
                    strainOptions=strainOptions,
                    pHLimit=pHLim,
                    pHVal=pHval,
                    networkAnalysis=TRUE
                                         #checkingOptions=list(stoiTol=0.1)
                )
                
                test.flag <<- testCase #change flag
            }
        }
            
        }else{#function has not been run before

            out <<- microPopModel(
                microbeNames=microbeNames,
                times=times,
                resourceSysInfo=resInfo,
                microbeSysInfo=micInfo,
                plotOptions=list(plotFig=FALSE),
                numStrains=numStrains,
                strainOptions=strainOptions,
                                        #checkingOptions=list(stoiTol=0.1)
                pHLimit=pHLim,
                pHVal=pHval,
                networkAnalysis=TRUE
           )
            
            test.flag <<- testCase #set flag
            
        }


    resourceNames <<- out$parms$resourceNames
    microbeNames <<- out$parms$microbeNames
    numPaths <<- out$parms$numPaths
    stoichiom <<- out$parms$Pmats$stoichiom
    Rtype <<- out$parms$Pmats$Rtype
    halfSat <<- out$parms$Pmats$halfSat
    maxGrowthRate <<- out$parms$Pmats$maxGrowthRate
    allStrainNames <<- out$parms$allStrainNames
    keyRes <<- out$parms$keyRes
    yield <<- out$parms$Pmats$yield
    parms <<- out$parms
    numStrains <<- out$parms$numStrains
}
