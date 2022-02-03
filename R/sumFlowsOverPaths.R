#' sumFlowsOverPaths
#'
#' sum flows over links between the same nodes
#' i.e. if the link has more than one metabolic path
#'
#' @param links data frame or matrix of links
#' @return matrix of links 
#' @export
#' 
sumFlowsOverPaths=function(links){
    
    nlinks=matrix(NA,ncol=ncol(links),nrow=nrow(links))
    colnames(nlinks)=colnames(links)
    
    ct=1
    for (i in 1:nrow(links)){

        start=links[i,'from']

        if (start!='Counted'){
            end=links[i,'to']

            ii=links[,'from']==start & links[,'to']==end
         
            nlinks[ct,]=links[i,]
            nlinks[ct,'weight']=sum(as.numeric(links[ii,'weight']))
            nlinks[ct,'path']='summed'

            links[ii,'from']='Counted'
        }
        
        ct=ct+1

    }

    sum.links=nlinks[!is.na(nlinks[,'from']),]

    return(sum.links)
}


    

    
