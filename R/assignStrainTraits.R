#' Internal function to assign stochastic strain traits
#' 
#' Produces a random distribution of trait values where the mean is groupVal
#' and the range is determined by strainOptions$percentTraitRange (if not pHtrait)
#' and by maxPHshift if it is the pHtrait (see strainOptions)
#' @param numStrains Integer. Number of strains per group
#' @param groupVal Scalar. Group parameter value i.e. the mean parameter value 
#' @param strainOptions list from microPopModel inputs. 
#' Contains 'distribution' i.e. the shape of the distribution ('normal' or 'uniform').
#' If it is not for a pH trait and the distribution is
#' 'normal' then its std dev is groupVal*percentRange/200, if distribution is
#' 'uniform' then its range is groupVal*(1 +/- percentRange/100).
#' For a pH trait, 'maxPHshift' is the max shift in pH units and
#' 'normal' has std dev = maxPHshift/2, and 'uniform' distribution has
#' range groupVal +/- maxPHshift;
#' @param parName Name of parameter. This is only used to help with error catching
#' @param pHtrait TRUE/FALSE whether or not trait is the pH trait.
#' @param gname Microbe name (for indexing strainOptions$percentTraitRange)
#' @importFrom stats rnorm runif
#' @return vector of values for each strain for one parameter
#' @export
assignStrainTraits = function(numStrains, groupVal, strainOptions,
    parName = "unspecified param", pHtrait = FALSE, gname='None') {

    
    if (!is.na(groupVal)) {


        #input checking
        if (strainOptions$distribution != "normal" & strainOptions$distribution != 
            "uniform") {
            stop(paste("MICROPOP ERROR: please specify the type of distribution (uniform or normal) from which to choose strain traits from for", parName))
        }
        if (length(strainOptions$percentTraitRange)>1){
            if (gname=='None'){
                stop(paste("MICROPOP ERROR: strainOptions$percentTraitRange has multiple numbers, please specify which group to use"))
            }
        }
        if (length(strainOptions$maxPHshift)>1){
            if (gname=='None'){
                stop(paste("MICROPOP ERROR: strainOptions$maxPHshift has multiple numbers, please specify which group to use"))
            }
        }
         
     #---------------------------------------

        if (pHtrait) {
            if (length(strainOptions$maxPHshift)>1){
                Range = strainOptions$maxPHshift[gname]
            }else{
                Range = strainOptions$maxPHshift
            }
        } else {
            if (length(strainOptions$percentTraitRange)>1){
                Range = groupVal*strainOptions$percentTraitRange[gname]/100
            }else{
                Range = groupVal*strainOptions$percentTraitRange/100
            }
        }


        if (strainOptions$distribution == "uniform") {

            
            vals = runif(numStrains, min = groupVal-Range, max = groupVal+Range)

            if (any(vals < 0)) {
                stop(paste("MICROPOP ERROR: The values of the mean and percent Trait Range for choosing strain trait values for", 
                  parName, "will result in values that are less than zero and therefore not valid"))
            }
        }

        
        if (strainOptions$distribution == "normal") {
            # draw parameter values from a distribution where 2*sigma=Range
            vals = rnorm(numStrains, mean = groupVal, sd = Range/2)
            if (any(vals < 0)) {
                stop(paste("MICROPOP ERROR: The values of the mean and std deviation for choosing strain trait values for", 
                  parName, "will result in values that are less than zero and therefore not valid"))
            }
        }
        
        
        if (any(numStrains > 1) & (is.na(Range) | Range == 0)) {
            print(paste("MICROPOP WARNING: There are", numStrains, "but they will all have the same value for", 
                        parName, "unless percentTraitRange is specified (and is above zero)"))
        }
        
        
        
    } else {
        vals = NA * seq(1, numStrains)
    }
    
    return(vals)
}
