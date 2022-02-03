#' Generic plotting of microbes over time
#' @param out output from microPopModel() 
#' @param yLabel String for y axis label. Default is 'Concentration'
#' @param xLabel String for x axis label. Default is 'Time'
#' @param sumOverStrains Logical. Default=TRUE
#' @param legendPosition String. Position of legend in microbe plot, default is 'topleft'
#' @param cex.title Scaling for title text
#' @param cex.ax Scaling for axes text (labels and ticklabels)
#' @param cex.legend  Scaling for legend text
#' @return Nothing just generates a plot
#' @importFrom graphics legend lines plot
#' @export
plotMicrobes = function(out,
    sumOverStrains=TRUE,
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
    microbeNames=out$parms$microbeNames

    numMFG = length(microbeNames)
    if (length(numStrains)==1){
        numM = numStrains * numMFG
    }else{
        numM =sum(numStrains)
    }
    
    time = soln[, 1]
    
    if (numM == 1) {
        # one strain and one group (need to keep as a matrix)
        X = as.matrix(soln[, 2:(numM + 1)])
        colnames(X) = colnames(soln)[2]
    } else {
        X = soln[, 2:(numM + 1)]
    }

    
    if (any(numStrains > 1) & sumOverStrains) {
        gmat = matrix(NA, nrow = length(time), ncol = numMFG)
        for (g in 1:numMFG) {
            gname=microbeNames[g]
            ssum=rep(0,length(time))
            for (strain in colnames(X)){
                if (!is.na(pmatch(gname,strain))){
                    ssum = ssum + X[,strain]
                }
            }
            gmat[, g] = ssum
        }
        Xmax = max(gmat, na.rm = TRUE)
    } else {
        Xmax = max(X, na.rm = TRUE)
    }
    

    cols = rainbow(numMFG)  #different groups have different colours
    plot(range(time), c(min(0, min(X, na.rm = TRUE)), 1.1 * Xmax), xlab = xLabel, 
         main = "Microbes", ylab = yLabel,
         cex.lab = cex.ax, cex.axis = cex.ax, cex.main = cex.title, 
         type = "n")
    for (g in 1:numMFG) {
        gname=microbeNames[g]
        if (any(numStrains > 1) & sumOverStrains) {
            lines(time, gmat[, g], lwd = 2, col = cols[g])
        } else {
            for (strain in colnames(X)){
                if (!is.na(pmatch(gname,strain))){
                    lines(time,X[,strain],lwd=2,col=cols[g])
                }
            }
        }
    }
    legend(legendPosition, bg = "transparent", legend = microbeNames, col = cols, lty = 1, 
           lwd = 2,cex=cex.legend)
}
