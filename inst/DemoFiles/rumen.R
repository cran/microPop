#rumen model
#3 MFGs: Xsu, Xaa, Xh2
#No dependence on pH

microbeNames=c('Xsu','Xh2','Xaa')
polymer.names=c('Zndf','Znsc','Zpro')

#Need to add polymers as extra resources c('Zndf','Znsc','Zpro') to a microbe dataframe so that they become state variables
Xsu[['Zndf']]=c('X',rep(NA,6))
Xsu[['Znsc']]=c('X',rep(NA,6))
Xsu[['Zpro']]=c('X',rep(NA,6))

time.hours=seq(0,3,1/(24*4))*24

#parameters
kd = 8.3333e-4 #death rate of microbes
f.X=c(0.2,0.55); names(f.X)=c('Znsc','Zpro') #fraction of dead cells that turn to polymers
khyd=c(0.05,0.20,0.22); names(khyd)=c('Zndf','Znsc','Zpro')#hydrolysis of polymers

#input files
resourceDF=createDF(paste(path.package('microPop'),'/inst/DemoFiles/resourceSysInfoRumen.csv',sep=''))
microbeDF=createDF(paste(path.package('microPop'),'/inst/DemoFiles/microbeSysInfoRumen.csv',sep=''))
Xh2['maxGrowthRate','Sh2']=0.6 #increase growth rate of methanogens

myRateFuncs=rateFuncsDefault

myRateFuncs$removalRateFunc=function(varName,varValue,stateVarValues,time,washOut,parms){
    death=0;    hydrolysis=0       
    if (varValue<=0){
        v=0
    }else{
        if (varName%in%polymer.names){hydrolysis=khyd[varName]}#hydrolysis of polymers
        if (getGroupName(varName,microbeNames)%in%microbeNames){death=kd}#death of microbes
        v=(washOut[varName] + hydrolysis + death)*varValue 
    }
    return(v)
}

myRateFuncs$entryRateFunc=function(varName,varValue,stateVarValues,time,inflowRate,parms){

    #entry rate from outside the system----------------------------------------------
    gname=getGroupName(varName,microbeNames)
    if (varName%in%parms$microbeNames){
        v.in=inflowRate[gname]/parms$numStrains
    }else if (varName%in%polymer.names){
        v.in=inflowRate[varName]*max(sin(pi*time/6),0)
    }else{
        v.in=inflowRate[varName]
    }

    #entry rate from inside the system-----------------------------------------------
    hydrolysis=0;      input.from.dead.cells=0
    
    #entry rate of sugar from hydrolysed polymers
    if (varName=='Ssu'){
        hydrolysis=sum(khyd[c('Zndf','Znsc')]*stateVarValues[c('Zndf','Znsc')])}#g/L/h

    #entry rate of amino acids from hydrolysed proteins
    if (varName=='Saa'){
        hydrolysis=khyd['Zpro']*stateVarValues['Zpro']}#g/L/h
  
    #dead microbial cells become Znsc and Zpro
    if (varName=='Znsc' | varName=='Zpro'){
        input.from.dead.cells=f.X[varName]*kd*sum(stateVarValues[parms$microbeNames])} #g/L/h

    v=v.in+hydrolysis+input.from.dead.cells
    return(v)
}

out=microPopModel(
    microbeNames=microbeNames,
    times=time.hours,
    rateFuncs=myRateFuncs,
    resourceSysInfo=resourceDF,
    microbeSysInfo=microbeDF,
  checkingOptions=list(balanceTol=1e-2,reBalanceStoichiom=FALSE,checkMassConv=FALSE,checkStoichiomBalance=FALSE),
    plotOptions=list(yLabel='concentration (g/l)',xLabel='time (h)',plotFig=TRUE,sumOverStrains=FALSE,saveFig=FALSE,figType='eps',figName='Rumen')
)

#plot methane concentration
dev.new()
par(mar=c(5,5,5,2))
plot(time.hours,out$solution[,'Sch4'],type='l',xlab='Time (hours)',ylab='Methane concentration (g/l)',cex.lab=1.5,cex.axis=1.5,lwd=2,main='Methane',cex.main=1.5)
#dev.copy2eps(file='RumenMethane.eps')
