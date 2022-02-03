
#' getVNPlotObject
#' 
#' uses visNetwork to produce an interactive network plot based on the links and edges dataframes
#'
#' @param nodes data frame or a list with nodes information. Needs at least column "id". See visNetwork::visNodes
#' @param edges data frame or a list with edges information. Needs at least columns "from" and "to". See visNetwork::visEdges
#' @param addLegend Logical. If true adds a legend to plot. Default is FALSE
#' @param addExport Logical. If true adds button to export fig from html plot
#' @param figType Type of export. One of "png" (default), "jpeg" or "pdf". Puts a button on the html plot
#' @param mainTitle Optional list containing "text" (string for plot title) and "style" (e.g. 'font-family:Times','font-family:Arial' etc).
#' @param subTitle Optional list containing "text" (string for plot subtitle) and "style" (e.g. 'font-family:Times','font-family:Arial' etc)
#' @param layoutSeed : NA. Random seed for the layout of the plot. To get identical plots set this to a number
#' @param scaleNodes Logical. If true the node sizes differ with concentration (in moles for resources and mass or concentration for microbes)
#' @param scaleEdges Logical. If true the edge sizes differ with the amount of moles flowing through them
#' @param microbeCol String for microbe node colour. Default is 'orange'
#' @param resourceCol String for resource node colour. Default is 'lightBlue'
#' @param productionCol String for production edge colour. Default is 'darkGrey'
#' @param uptakeCol String for uptake edge colour. Default is 'magenta'
#' @param figWidth numeric value to control size of plotting window. Default is 700
#' @param figHeight numeric value to control size of plotting window. Default is 700

#' @import visNetwork 
#' @export
#' @return a visNetwork object that can be shown using print() function.

getVNPlotObject=function(nodes,edges,
                         addLegend=FALSE,
                         addExport=TRUE,
                         figType='png',
                         mainTitle=NULL,
                         subTitle=NULL,
                         layoutSeed=NA,
                         scaleNodes=FALSE,
                         scaleEdges=FALSE,
                         microbeCol='gold',
                         resourceCol='lightblue',
                         productionCol='magenta',
                         uptakeCol='darkgrey',
                         figWidth=700,
                         figHeight=700){

    if (is.null(nodes) | is.null(edges)){
        stop("MICROPOP ERROR: Must have 'nodes' and 'edges' data frames input arguments for getVNPlotObject()")
    }

    if (scaleNodes){
        #scale microbes and resources nodes differently
        ii=seq(1,nrow(nodes))[nodes[,'group']=='microbe']
        microbe.vals=as.numeric(nodes$value[ii])/max(as.numeric(nodes$value[ii]),na.rm=TRUE)

        iii=seq(1,nrow(nodes))[nodes[,'group']=='resource']
        res.vals=as.numeric(nodes$value[iii])/max(as.numeric(nodes$value[iii]),na.rm=TRUE)

        nodes$value[ii]=microbe.vals
        nodes$value[iii]=res.vals
        scaling=list(min=min(nodes[,'value']),max=max(nodes[,'value']))
        
    }else{
        scaling=NULL
        nodes[,'value']=1
    }

    if (scaleEdges){
        widths=as.numeric(edges[,'width'])
        w=(widths-min(widths))/diff(range(widths))
        edges[,'width']=1+7*w
    }else{
        edges[,'width']=2
    }
    
    node.cols=rep(NULL,nrow(nodes))
    node.cols[nodes[,'group']=='microbe']=microbeCol
    node.cols[nodes[,'group']=='resource']=resourceCol

    nodes$color=node.cols
    
    edge.cols=rep(NULL,nrow(edges))
    edge.cols[edges[,'group']=='uptake']=uptakeCol
    edge.cols[edges[,'group']=='production']=productionCol

    edges$color=edge.cols

    vv=visNetwork(nodes, edges,                              
                  width=figWidth,height=figHeight,
                  main=mainTitle,submain=subTitle,
                  scaling=scaling
                  )
  
    if (addExport){
        vv = vv %>% visExport(type=figType)
    }
    
    if (addLegend){

        ledges <- data.frame(color = c(uptakeCol,productionCol),
                             label = c("uptake", "production"), 
                             arrows = c("to", "from"),
                             width=2,stringsAsFactors=FALSE)
        
        vv=vv%>%
            visGroups(groupname = "microbe", color = microbeCol) %>%
            visGroups(groupname = "resource", color = resourceCol) %>%
            visLegend(addEdges=ledges,useGroups = TRUE)
    }
    
    if (!is.na(layoutSeed)){
        vv = vv %>% visLayout(randomSeed=layoutSeed)
    }
    
    return(vv)
}
