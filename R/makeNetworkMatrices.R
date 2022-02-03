#' makeNetworkMatrices
#'
#' make links and nodes matrices for use in network plotting software 
#'
#' @param chosen.time the time you want to plot
#' @param out the output from microPopModel()
#' @param convertToMoles Logical. Default is TRUE
#' @param sumOverStrains Logical. Default is TRUE
#' @export
#' 
makeNetworkMatrices = function(chosen.time,out,convertToMoles=TRUE,sumOverStrains=TRUE) {

    #requires functions: reshapeFlowMat, convertFlowsToMoles, convertStatesToMoles
    if (chosen.time<min(out$solution[,'time']) | chosen.time>max(out$solution[,'time'])){
        stop('MICROPOP ERROR: chosen.time must be within the simulation time')
    }
    
    #find time step closest to chosen.time
    time.step=stats::approx(out$solution[,'time'],seq(1,nrow(out$solution)),chosen.time)$y
 
    uptakeList=reshapeFlowMat(time.step,'uptake',out)
    productionList=reshapeFlowMat(time.step,'production',out)
    
    groupNames=out$parms$microbeNames
    allStrainNames=out$parms$allStrainNames
    
    resourceNames=out$parms$resourceNames
    rNodes=out$parms$resourceNames

    concentration.orig=out$solution[time.step,]
    
    if (sumOverStrains & length(groupNames)!=length(allStrainNames)){ #i.e. multiple strains in group

        uptakeList.new=sumFlowOverStrains(uptakeList,allStrainNames,groupNames)
        productionList.new=sumFlowOverStrains(productionList,allStrainNames,groupNames)
        
        concentration=sumConcOverStrains(concentration.orig,allStrainNames,groupNames,resourceNames)
        
    }else{
        
        uptakeList.new=uptakeList
        productionList.new=productionList
        concentration=concentration.orig
    }
    
        #print(productionList.new)
    
    if (convertToMoles){
        
        #convertFlowsToMoles
        uptake=convertFlowsToMoles(groupNames,
                                   uptakeList.new,out$parms$molarMass)
        production=convertFlowsToMoles(groupNames,
                                       productionList.new,out$parms$molarMass)

        #convertStatesToMoles
        nodeWeights=convertStatesToMoles(concentration,out$parms$molarMass)
        
    }else{

        uptake=uptakeList.new
        production=productionList.new
        nodeWeights=concentration

    }

    
    Links = matrix(NA, nrow = 1, ncol = 6)
    colnames(Links) = c('from', 'to', 'path','weight','type','microbe')
    
    if (sumOverStrains){
        microbeNodeNames=groupNames
    }else{
        microbeNodeNames=allStrainNames
    }
    
    for (microbe in microbeNodeNames){ 
        
        uptake.p=uptake[[microbe]]
        production.p=production[[microbe]]

        
        for (path in rownames(uptake.p)){
                
            for (sub in colnames(uptake.p)){

                if (!is.na(uptake.p[path,sub])){
                    if (uptake.p[path,sub]>0){
                        
                        Links=rbind(Links,c(sub,microbe,path,uptake.p[path,sub],'uptake',microbe))
                    }
                }
            }
            
            
            for (prod in colnames(production.p)){
                
                if (!is.na(production.p[path,prod])){
                    if (production.p[path,prod]>0){
                        
                        Links=rbind(Links,c(microbe,prod,path,production.p[path,prod],'production',microbe))
                    }
                }
            }
                
        }#path
    }#microbe
    
    links = Links[2:nrow(Links), ]

    nodes = cbind('id'=c(rNodes,microbeNodeNames),'type'=c(rep('resource',length(rNodes)),rep('microbe',length(microbeNodeNames))),'concentration'=nodeWeights[c(rNodes,microbeNodeNames)])

    return(list(links=links,nodes=nodes))
    
}
