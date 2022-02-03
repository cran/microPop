#'
#' networkDFfromMPinput
#'
#' make node and edge data frames to use in visNetwork from microPop microbial data frames
#'
#' @param microbeNames vector of strings of the names of the microbial data frames you want to plot. These can be intrinsic data frames or loaded in by user.
#' 
#' @export
#' @return a list containing the edges and nodes


networkDFfromMPinput=function(microbeNames){

    mfgList = NULL
    
    for (g in 1:length(microbeNames)) {
      
      gname = microbeNames[g]

      group = get(gname)
      
         # get names of substrates the microbe uses:

        subs=NULL
        prods=NULL

            txt='Rtype'
            i=1

            while (txt%in%rownames(group)){

                newsubs=colnames(group)[
                    group[txt, ] == 'S' |
                    group[txt, ] == 'Se' |
                    group[txt, ] == 'Sw' |
                    group[txt, ] == 'Sb' ]
                
                subs=c(subs,newsubs)
                
                prods=c(prods,colnames(group)[group[txt, ] == 'P' ])

                i=i+1
                txt=paste('Rtype',i,sep='.')

            }
          
        mfgList[[gname]][['substrate']] =unique(subs)
        mfgList[[gname]][['product']] =unique(prods)
    }
    
    rNodes = getAllResources(microbeNames)
    XNodes = microbeNames

    links=rep(NULL,3)
     for (x in XNodes) {
      for (r in rNodes) {
        if (r %in% mfgList[[x]]$'substrate') {
            links=rbind(links,c(r, x, 'uptake'))
        }
        if (r %in% mfgList[[x]]$'product') {
            links=rbind(links,c(x, r, 'production'))
        }
      }
    }
    
    colnames(links) = c('from', 'to', 'type')
    
    edges <- data.frame(from = links[,'from'], to = links[,'to'],arrows='to',
                        group=links[,'type'],
                        width=1,
                        stringsAsFactors=FALSE)
    
 
    Nodes = c(rNodes, XNodes)
    Node.type = c(rep('resource', length(rNodes)), rep('microbe', length(XNodes)))
    nodes = data.frame(id=Nodes,
                       label=Nodes,
                       group=Node.type,
                       shadow=TRUE,
                       value=1,
                       stringsAsFactors=FALSE)

    return(list(nodes=nodes,edges=edges))
    
}



                         
