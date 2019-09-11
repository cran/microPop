#' Makes vector of unique resource names
#' @param microbeNames Vector of strings which contains the names of the microbial groups in the system e.g. c('Bacteroides','Acetogens')
#' @return vector of resource names
#' @export
#' 
getAllResources = function(microbeNames) {
    
    # read in all the microbe files and get the full list of resources needed
    
    allResources = NULL

    for (gname in microbeNames) {
        data = get(gname)
        nres = colnames(data)
        for (i in 1:length(nres)){
            if (!is.na(data['Rtype',i]) & nres[i] != "units"  & nres[i] != "Units" &
                nres[i] != "Biomass" & nres[i] != "biomass" ){
                allResources=c(allResources,nres[i])
            }
        }
    }
    
    resNames = unique(allResources[!is.na(allResources)])
    
    return(resNames)
    
}
