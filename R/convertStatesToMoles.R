#' convertStatesToMoles
#'
#' convert network nodes from mass to moles for resources (microbes remain as mass)
#'
#' @param nodeMass is the value of each node in the network (named vector)
#' @param MolarMass is a named vector containing the molar mass for each resource e.g. out$parms$molarMass
#' @export
#' 
convertStatesToMoles=function(nodeMass,MolarMass){

    nodeWeights=nodeMass
    
    #convert resources to moles
    for (s in names(nodeWeights)){
        if (s%in%names(MolarMass)){
            nodeWeights[s]=nodeMass[s]/MolarMass[s]
        }
    }
    
    return(nodeWeights)

}#----------------------------------------------------------------------------------------------------------------------------------------
