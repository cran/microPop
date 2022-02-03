#' get stochastically generated pH corners for each strain
#' 
#' Returns the values of the pH values of the limit function i.e. where the limit is c(0,1,1,0)
#' Reads these in from the microbe group dataframes
#' 
#' @param microbeNames (vector of strings). Names of microbes in the system
#' @param allStrainNames (vector of strings)
#' @param numStrains Integer or named vector of integers
#' @param pHcorners vector of 4 scalars definining the pH lim func
#' @param pHLimit (logical) Is microbial growth affected by pH?
#' @param strainOptions list from microPopModel inputs
#' @param oneStrainRandomParams logical from microPopModel inputs
#'
#' @return (matrix) values of the pH values of the limit function i.e.
#' where the limit is c(0,1,1,0) for each strain
#' @export


getStrainPHcorners = function(microbeNames, allStrainNames, numStrains, pHcorners, 
    pHLimit, strainOptions,oneStrainRandomParams) {

#    mat = matrix(NA, ncol = 4, nrow = length(allStrainNames),
#        dimnames = list(allStrainNames))

    mat=NULL
    
    if (pHLimit){
        for (g in 1:length(microbeNames)) {

            #check if numStrains has a value per group
            if (length(numStrains)==1){
                Ls=numStrains
            }else{
                Ls=numStrains[microbeNames[g]]
            }

            #print(strainOptions)
            shifts=rep(0,Ls)
            
            if (Ls>1 | oneStrainRandomParams){
                
                if ("pHtrait" %in% strainOptions$randomParams) {
                    shifts = assignStrainTraits(Ls, 1, strainOptions,
                                                parName = "pHtrait", 
                                                pHtrait = TRUE, microbeNames[g]) - 1
                }
            }
            
            for (i in 1:Ls) {
                mat=rbind(mat,pHcorners[microbeNames[g], ]+shifts[i])
            }
        }

        rownames(mat)=allStrainNames
        
    }

       
    return(mat)
}


