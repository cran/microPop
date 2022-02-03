#' reshapeFlowMat
#'
#'reshapes the flow matrices out$flow.uptake or out$flow.production into a list
#' elements of the list are the microbeNames and then there is a matrix [path,res]
#' 
#' @param time.step is the index of the chosen time
#' @param flow.direction is either 'uptake' or 'production'
#' @param out is the output from microPopModel with networkAnalysis=TRUE
#' @return a list with microbeNames as elements and a matrix of [path,resource] showing the chosen flow direction (eg. uptake or production). Note theses flows have not been converted to moles.
#' 
#' @export
#' 
reshapeFlowMat=function(time.step,flow.direction,out){

    if (flow.direction=='uptake' | flow.direction=='Uptake' ){
        if ('flow.uptake'%in%names(out)){
            flow.mat=out$flow.uptake
        }else{
            stop('MICROPOP ERROR: the flow output from microPop is missing. Please add networkAnalysis=TRUE to microPopModel input arguments')
        }
    }else if (flow.direction=='production' | flow.direction=='Production' ){
        if ('flow.production'%in%names(out)){
            flow.mat=out$flow.production
        }else{
            stop('MICROPOP ERROR: the flow output from microPop is missing. Please add networkAnalysis=TRUE to microPopModel input arguments')
        }
    }
    

    numR=length(out$parms$resourceNames) #number of resources
    numS=length(out$parms$allStrainNames) #number of strains
    numP=max(out$parms$numPaths)  #max number of paths
    allStrainNames=out$parms$allStrainNames

    flow.row=flow.mat[time.step,]
    flowList=list(allStrainNames)

    for (bac in allStrainNames){

        new.mat=matrix(NA,ncol=numR,nrow=numP)
        colnames(new.mat)=out$parms$resourceNames
        rownames(new.mat)=paste('path',seq(1,numP),sep='')

        for (res in out$parms$resourceNames){
            strng=paste(bac,'.',res,'.path',1:max(out$parms$numPaths),sep='')
            new.mat[,res]=flow.row[strng]
        }

        flowList[[bac]]=new.mat

    }

    return(flowList)
    
}
