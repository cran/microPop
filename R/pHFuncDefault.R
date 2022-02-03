#' pH Function
#'
#' Return the value of pH in pH units
#'
#' @aliases pHFunc
#' @param time (scalar). The current time point in the ODE solver.
#' @param parms List which contains all information required by the ODE solver
#' @param stateVarValues State vector (resources and microbes) (with names)
#' @return (scalar) pH at the given time
#' @export
pHFuncDefault <- function(time, parms, stateVarValues=NULL) {
    pH = parms$pHVal
    
    return(pH)
}
