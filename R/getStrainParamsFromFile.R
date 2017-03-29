#' get strain parameter values from a csv file
#' 
#' @param Pmats List of parameter matrices
#' @param strainPHcorners Matrix of pH corners for each strain
#' @param strainOptions List which is input to microPopModel
#' @return (list) - first entry is new version of Pmats, second is new version of strainPHcorners

getStrainParamsFromFile=function(Pmats,strainPHcorners,strainOptions){

    nPmats=Pmats
    nCorners=strainPHcorners

    filename=strainOptions$paramFileName
    data=read.csv(filename,header=TRUE,stringsAsFactors=FALSE,row.names=NULL)

    for (line in 1:length(data[,1])){
        strainName=trimws(data[line,'strainName'],which='both')
        path=as.numeric(trimws(data[line,'path'],which='both'))
        resName=trimws(data[line,'resource'],which='both')
        parVal=as.numeric(trimws(data[line,'paramVal'],which='both'))
        parName=trimws(data[line,'paramName'],which='both')
        if (parName=='pHcorners'){
            nCorners[strainName,]=as.numeric(data[line,3:6])
        }else{
            nPmats[[parName]][[strainName]][path,resName]=parVal
        }
    }

    return(list(nPmats,nCorners))
}

