#'
#' networkDFfromMPoutput
#'
#' make node and edge data frames from microPop output to use in visNetwork
#'
#' @param chosen.time the time you want to plot
#' @param MPoutput the output from microPopModel()
#' @param groupNames Default is NULL which plots all the microbes. To plot a subset of all the groups, specify a vector of strings of the names of the groups you want to plot.
#' @param sumOverPaths  Logical. Default is TRUE which sums flows between the same nodes even if they are on different metabolic paths
#' @param sumOverStrains Logical. Default is TRUE which means the strains are put into their functional group nodes and the flow are summed. When it is FALSE, each strain will have its own node.
#' @param convertToMoles Logical. Default is TRUE

#' @export
#' @return a list containing the edges and nodes


networkDFfromMPoutput=function(chosen.time,
                               MPoutput,
                               groupNames=NULL,
                               sumOverPaths=TRUE,
                               sumOverStrains=TRUE,
                               convertToMoles=TRUE
                               ){

    #requires function: sumFlowsOverPaths and makeNetworkMatrices which in turn requires:
    #reshapeFlowMat, convertFlowsToMoles and convertStatesToMoles
    
    networkObj=makeNetworkMatrices(chosen.time,MPoutput,convertToMoles,sumOverStrains)
    linksMP=networkObj$links
    nodesMP=networkObj$nodes

    if (sumOverPaths){
        #sum flows which have multiple metabolic paths
        linksMP=sumFlowsOverPaths(linksMP)
    }

    #subset of groups
    if (!is.null(groupNames)){
        keepE=NULL
        for (i in 1:nrow(linksMP)){
            if (linksMP[i,'from']%in%groupNames | linksMP[i,'to']%in%groupNames){
                keepE=c(keepE,i)
            }
        }
    }else{
        keepE=seq(1,nrow(linksMP))
    }

    
    edges=data.frame(
        from=linksMP[keepE,'from'],
        to=linksMP[keepE,'to'],
        width=linksMP[keepE,'weight'],
        arrows='to',
        group=linksMP[keepE,'type'],
        stringsAsFactors=FALSE
    )


    #subset of groups
    if (!is.null(groupNames)){
        keepN=NULL
        for (i in 1:nrow(nodesMP)){
            if (nodesMP[i,'id']%in%edges[,'from'] | nodesMP[i,'id']%in%edges[,'to']){ 
                keepN=c(keepN,i)
            }
        }
    }else{
        keepN=seq(1,nrow(nodesMP))
    }

    
    nodes=data.frame(
        id=nodesMP[keepN,'id'],
        label=nodesMP[keepN,'id'],
        shadow=TRUE,
        group=nodesMP[keepN,'type'],
        value=as.numeric(nodesMP[keepN,'concentration']),
        stringsAsFactors=FALSE
    )


    return(list(edges=edges,nodes=nodes))

}
