#' Checks whether the solution generated by the ODE solver contains negative values
#' 
#' @param soln Matrix from ode solver out$solution
#' @param tol tolerance
#' @export

checkSolution = function(soln, tol = -0.1) {

    stateVarNames=colnames(soln)
    for (v in 2:length(stateVarNames)) {
        var = stateVarNames[v]
        if (any(soln[, var] < tol,na.rm=TRUE)) {
            warning(paste("MICROPOP WARNING: There are negative values for", var))
        }
        
    }
}
