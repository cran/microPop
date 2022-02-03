#' convertFlowsToMoles
#'
#' convert network flows from mass to moles 
#'
#' @param allStrainNames is a vector containing the names of the microbes (strings) 
#' @param flow is the list output from reshapeFlowMat()
#' @param molarMass is a named vector containing the molar mass for each resource e.g. out$parms$molarMass
#' @export
#' 
convertFlowsToMoles=function(allStrainNames,flow,molarMass){

    flow.moles=flow
    
    for (g in allStrainNames){

        mat=flow[[g]]
        mat.moles=mat

        for (s in colnames(mat)){
            if (s%in%names(molarMass)){
                mat.moles[,s]=mat[,s]/molarMass[s]
            }else{
                stop(paste('Molar mass for',s,'is not in',molarMass))
            }
        }
        flow.moles[[g]]=mat.moles

    }

    return(flow.moles)
    
}
