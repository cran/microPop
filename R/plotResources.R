#' Generic plotting of resources over time
#' @param out output from microPopModel()
#' @param yLabel String for y axis label. Default is 'Concentration'
#' @param xLabel String for x axis label. Default is 'Time'
#' @param legendPosition String. Position of legend in resource plot, default is 'topleft'
#' @param cex.title Scaling for title text
#' @param cex.ax Scaling for axes text (labels and ticklabels)
#' @param cex.legend Scaling for legend text
 #' @return Nothing just generates a plot
#' @importFrom graphics legend lines par plot
#' @export
plotResources = function(out,
    yLabel='Concentration',
    xLabel='Time',
    legendPosition="topleft",
    cex.title=1,
    cex.ax=1,
    cex.legend=1
                         ){

    soln=out$solution
    numR=length(out$parms$resourceNames)
    numStrains=out$parms$numStrains
    
    numMFG = length(out$parms$microbeNames)
    if (length(numStrains)==1){
        numM = numStrains * numMFG
    }else{
        numM =sum(numStrains)
    }
    
    time = soln[, 1]
    
    R = soln[, (numM + 2):(numM + numR + 1), drop = FALSE]
    
    cols = rainbow(numR)
    plot(range(time), c(min(0, min(R, na.rm = TRUE)), 1.1 * max(R, na.rm = TRUE)), 
        xlab = xLabel, main = "Resources", ylab = yLabel,
         cex.lab = cex.ax, cex.axis = cex.ax, 
        type = "n", cex.main = cex.title)
    for (i in 1:numR) {
        lines(time, R[, i], lwd = 2, col = cols[i])
    }
    legend(legendPosition, bg = "transparent", colnames(R), col = cols, lty = 1, lwd = 2,cex=cex.legend)
}
