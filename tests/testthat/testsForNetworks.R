#tests for network functions

#Note loadTestDataFunc() is called in each test (so they are
#stand alone) but microPopModel() will only be run if it hasn't
#already been run before

source('loadTestDataFunc.R')


test_that('Test convertStatesToMoles.R',{

    #does not require a model run
    
    nm=c('one','two','three')
    nodeMass=c(1,2,3)
    names(nodeMass)=nm
    MolarMass=c(2,1,3)
    names(MolarMass)=nm
    v=convertStatesToMoles(nodeMass,MolarMass)

    res=c(0.5,2.0,1.0)
    names(res)=nm
    expect_equal(v,res)
    
})

test_that('Test convertFlowsToMoles.R',{

    #requires a model run

    for (case in 1:2){
        loadTestDataFunc(case)
        times=out$solution[,'time']
        flow=reshapeFlowMat(round(length(times)/2),
                            flow.direction='uptake',
                            out)
        
        flow.moles=convertFlowsToMoles(out$parms$allStrainNames,flow,out$parms$molarMass)
        strain1=out$parms$allStrainNames[1]
        
                                        #check that the molar flow is not identical to the original flow
        expect_false(identical(flow[[strain1]],flow.moles[[strain1]]))
        
                                        #check that molar flow is mass flow/molar mass
        expect_equal(flow[[strain1]][1,]/out$parms$molarMass[colnames(flow[[strain1]])],flow.moles[[strain1]][1,])
    }
    
    
})


test_that('Test reshapeFlowMat.R',{

    #requires a model run

    for (case in 1:2){
        loadTestDataFunc(case)
        times=out$solution[,'time']
        v=reshapeFlowMat(round(length(times)/2),
                         flow.direction='uptake',
                         out)
        
        expect_equal(v[[1]],out$parms$allStrainNames)
        
        strain1=out$parms$allStrainNames[1]
        group1=getGroupName(strain1,out$parms$microbeNames)
        
                                        #check all the resources for the group are included in the matrix
        expect_true(all(getAllResources(group1)%in%colnames(v[[strain1]])))
        
        strain2=out$parms$allStrainNames[length(out$parms$allStrainNames)]
        group2=getGroupName(strain2,out$parms$microbeNames)
        
                                        #check all the resources for the group are included in the matrix
        expect_true(all(getAllResources(group2)%in%colnames(v[[strain2]])))
    }
})


test_that('Test makeNetworkMatrices.R',{

    for (case in 1:2){
        loadTestDataFunc(case)
        times=out$solution[,'time']
        chosen.time=times[round(length(times)/2)]
        v=makeNetworkMatrices(chosen.time,out,convertToMoles=TRUE,sumOverStrains=TRUE)
        
        expect_true(all(unique(out$parms$microbeNames)%in%unique(v$links[,'to'])))

        expect_true(all(unique(v$links[,c('from','to')])%in%unique(c(out$parms$resourceNames,out$parms$microbeNames))))
        
        expect_true(all(unique(out$parms$microbeNames)%in%unique(v$links[,c('from','to')])))
 
        expect_true(sum(as.numeric(v$link[,'weight']))>=0)
        expect_true(sum(as.numeric(v$node[,'concentration']))>=0)
        expect_true(length(names(v))==2)

     }
    
})


test_that('Test sumFlowOverStrains.R',{

    mat1=0.1*rbind(path1=c(1,2,3),path2=c(4,5,6))
    mat2=0.2*rbind(path1=c(1,2,3),path2=c(4,5,6))
    mat3=0.3*rbind(path1=c(1,2,3),path2=c(4,5,6))
    mat4=0.4*rbind(path1=c(1,2,3),path2=c(4,5,6))
    mat5=0.5*rbind(path1=c(1,2,3),path2=c(4,5,6))
    mat6=0.6*rbind(path1=c(1,2,3),path2=c(4,5,6))

    allStrainNames=c(paste('A',1:3,sep='.'),paste('B',1:3,sep='.'))
    groupNames=c('A','B')
    
    flow=list(A.1=mat1, A.2=mat2, A.3=mat3,
              B.1=mat4, B.2=mat5, B.3=mat6)
              
    v=sumFlowOverStrains(flow,allStrainNames,groupNames)

    expect_equal(v$A,mat1+mat2+mat3)
    expect_equal(v$B,mat4+mat5+mat6)
    
})


test_that('Test sumConcOverStrains.R',{

    conc.orig=c(A.1=1,A.2=2,A.3=3,B.1=4,B.2=5,B.3=6,res1=1,res2=2)
    
    resourceNames=c('res1','res2')

    allStrainNames=c(paste('A',1:3,sep='.'),paste('B',1:3,sep='.'))
    groupNames=c('A','B')
    
    v=sumConcOverStrains(conc.orig,allStrainNames,groupNames,resourceNames)

    ans=c(res1=1,res2=2,A=6,B=15)
    expect_identical(v,ans)
     
})



test_that('Test sumFlowsOverPaths.R',{

    links1=cbind(from=c('A','A','A','C','C'),to=c('B','B','C','B','B'),weight=seq(1,5),path=c(1,2,1,1,2))

    v=sumFlowsOverPaths(links1)
    ii=(v[,'from']=='A' & v[,'to']=='B')
    iii=(v[,'from']=='A' & v[,'to']=='C')
    iiii=(v[,'from']=='C' & v[,'to']=='B')
    expect_equal(as.numeric(v[ii,'weight']),3)
    expect_equal(as.numeric(v[iii,'weight']),3)
    expect_equal(as.numeric(v[iiii,'weight']),9)

    expect_true(nrow(v)==3)

    
})

test_that('Test networkDFfromMPinput.R',{

     MFG=matrix(NA,ncol=4,nrow=6,dimnames=list(c('Rtype','halfSat','yield',
     'maxGrowthRate','stoichiom','keyResource'),c('H2','CO2','CH4','H2O')))
     MFG['Rtype',]=c('Se','Se','P','P')
     MFG['halfSat',c('H2','CO2')]=1e-6
     MFG['yield','H2']=0.2
     MFG['maxGrowthRate','H2']=2
     MFG['keyResource',1]='H2'
     MFG['stoichiom',]=c(4,1,1,2)
     Archea=data.frame(MFG,stringsAsFactors=FALSE)

     microbeNames='Archea'
     v=networkDFfromMPinput(microbeNames)

     expect_true(all(c('H2','CO2','CH4','H2O')%in%c(v$nodes[,'id'])))
     expect_true(all(c('H2','CO2')%in%v$edges[v$edges[,'to']==microbeNames,'from']))
     expect_true(all(c('CH4','H2O')%in%v$edges[v$edges[,'from']==microbeNames,'to']))


})



test_that('Test getVNPlotObject.R',{

    #this last test is very comprehensive
    #It is for a simple system with 1 microbe, 1 substrate and 1 product (case=3 in loadTestDataFunc)
    #We check microPop against the analytical solution and the deSolve solution
    
    Gmax=10
    K=0.01
    V=2
    Y=0.3
    SVin=10
    
    init=c(S=1,B=1,P=0)
    times=seq(0,10,0.001)
    
#call ODE solver from deSolve-------------------------------------------

    func=function(Time,y,parms){
        S=y[1]
        B=y[2]
        P=y[3]
        G=Gmax*S/(S+K)
        dB=B*G-V*B
        dS=SVin-B*G/Y-S*V
        dP=B*G*(1/Y-1)-P*V
        return(list(c(dS,dB,dP)))
    }

    outDS=ode(y=init,times=times,func=func,parms=NULL)
#Steady state - analytical solution-----------------------------------------
    Bss=Y*SVin/V
    Sss=0
    Pss=Bss*(1/Y-1)
    
#microPop solution---------------------------------------------------------------------------
    loadTestDataFunc(3)
    outMP=out
    
#flows
    S=outMP$solution[nrow(outMP$solution),'S']
    B=outMP$solution[nrow(outMP$solution),'B']
    S.to.B=(B*Gmax*S/(S+K))/Y
    B.to.P=B*(Gmax*S/(S+K))*(1/Y-1)
    
    DFs=networkDFfromMPoutput(chosen.time=max(times),MPoutput=outMP)
    
    vv=getVNPlotObject(DFs$nodes,DFs$edges,
                       addLegend=FALSE,
                       figWidth=700,
                       figHeight=700,
                       scaleNodes=TRUE,
                       scaleEdges=TRUE)
    #-----------------------------------------------------------------------------------------------------------------
    #TESTS
    #check microPop solution and analytic solution are the same
    expect_true(abs(sum(outMP$solution[nrow(outMP$solution),c('time','S','B','P')]-c(max(times),Sss,Bss,Pss)))<1e-5)
    #check microPop solution and deSolve solution are the same
    expect_true(abs(sum(outMP$solution[nrow(outMP$solution),c('time','S','B','P')]-outDS[nrow(outDS),]))<1e-10)

    #check nodes have correct concentration
    expect_true(abs(sum(outMP$solution[nrow(outMP$solution),c('S','B','P')]-DFs$nodes[c('S','B','P'),'value']))<1e-10)

     #Compare flows from microPop and analytical solution
    expect_true(S.to.B-as.numeric(DFs$edges[(DFs$edges[,'from']=='S' & DFs$edges[,'to']=='B'),'width'])<1e-5)
    expect_true(B.to.P-as.numeric(DFs$edges[(DFs$edges[,'from']=='B' & DFs$edges[,'to']=='P'),'width'])<1e-5)

    #check size ratios
    #resources are scaled separately to microbes
    expect_true(abs(vv$x$nodes['S','value']/vv$x$nodes['P','value']-outMP$solution[nrow(outMP$solution),'S']/outMP$solution[nrow(outMP$solution),'P'])<1e-5)

    #one of the resources must have size 1
    expect_true(any(vv$x$nodes[c('S','P'),'value']==1))
    #there is only one microbe so should have value 1
    expect_true(vv$x$nodes['B','value']==1)

})
