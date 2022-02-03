#' sumFlowOverStrains
#'
#' make links and nodes matrices for use in network plotting software 
#'
#' @param flowList is list containing the production or uptake flows (the output from reshapeFlowMat())
#' @param allStrainNames is a vector containing the names of the microbial strains (strings) 
#' @param groupNames is a vector containing the names of the microbial groups (strings)
#' @export
#' 
sumFlowOverStrains=function(flowList,allStrainNames,groupNames){

    flowList.new=list(groupNames)
    
    for (strain in allStrainNames){
        group=getGroupName(strain,groupNames) 
        if (strain==paste(group,'1',sep='.')){
            group.mat=flowList[[strain]]
        }else{
            group.mat=group.mat+flowList[[strain]]
        }
        flowList.new[[group]]=group.mat
    }

    return(flowList=flowList.new)

    
}
