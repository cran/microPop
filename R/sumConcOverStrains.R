#' sumConcOverStrains
#'
#' sum concentration of each strain into the group it is in
#'
#' @param concentration.orig the row of out$solution at the required time point
#' @param allStrainNames is a vector containing the names of the microbial strains (strings) 
#' @param groupNames is a vector containing the names of the microbial groups (strings)
#' @param resourceNames is a vector of strings containing the names of all the resources
#' 
#' @export
#' 
sumConcOverStrains=function(concentration.orig,allStrainNames,groupNames,resourceNames){

    concentration=rep(NA,length(groupNames)+length(resourceNames))
    names(concentration)=c(resourceNames,groupNames)
    concentration[resourceNames]=concentration.orig[resourceNames]
    
    for (strain in allStrainNames){
        group=getGroupName(strain,groupNames) 
        if (strain==paste(group,'1',sep='.')){
            group.conc=concentration.orig[strain]
        }else{
            group.conc=group.conc+concentration.orig[strain]
        }
        concentration[group]=group.conc
    }

    return(concentration)
    
}
