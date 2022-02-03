#' calculate the mean trait at the end of the model run
#' @param out Output from microPopModel()
#' @param trait.name can be 'halfSat','yield','maxGrowthRate' and 'pHtrait' or 'strainpHcorners'
#' @param gname name of group or microbe
#' @param resource.name String
#' @param path String
#' @export
#' 
meanTraitFunc = function(out, trait.name, gname, resource.name,path){
    
    if (length(out$parms$numStrains)>1){
        numStrains = out$parms$numStrains[gname] 
    }else{
        numStrains = out$parms$numStrains 
    }       
    
    if (numStrains<=1){ stop('Need more than one strain to look at trait changes')}

    
    strain.names = paste(gname, ".", seq(1, numStrains), sep = "")
    weighted.vals = NA*seq(1,numStrains); names(weighted.vals)=strain.names
        
    for (strain in strain.names) {
        if (trait.name == "strainPHcorners" | trait.name == "pHtrait") {
            trait.val = pHcentreOfMass(strain, gname, out$parms$pHLimFunc, out$parms)
        } else {
            trait.val = out$parms$Pmats[[trait.name]][[strain]][path, resource.name]
        }
        #print(trait.val)
        mass = out$solution[nrow(out$solution), strain]
        weighted.vals[strain] = trait.val * mass
    }
        
        # at each point in time find the sum(xi*Mi)/sum(Mi) over all strains i
    avTrait =sum(weighted.vals)/sum(out$solution[nrow(out$solution), strain.names])
        
    return(avTrait)
}



