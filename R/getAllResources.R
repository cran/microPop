#' Makes vector of unique resource names
#' @param microbeNames Vector of strings which contains the names of the microbial groups in the system e.g. c('Bacteroides','Acetogens')
#' @param gutModel Logical. TRUE if using with the microPopGut package
#' @param myPars list of extra parameters
#' @return vector of resource names
#' @export
#' 
getAllResources = function(microbeNames,gutModel=FALSE,myPars=NULL) {
    
    # read in all the microbe files and get the full list of resources needed
    
    allResources = NULL
    #ct = 0
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

    if (gutModel & !is.null(myPars)){
        if (!myPars$waterName%in%resNames){
            print(paste('adding',myPars$waterName,'to resourceNames'))
            resNames=c(resNames,myPars$waterName)
        }
    }
  
    return(resNames)
}
