#' Convert strain name to its group name
#' e.g. 'Bacteroides.1' becomes 'Bacteroides'
#' updated (Dec 2019) so that MFG names can contain dots 
#' @param xname a string (may be strain name or something else)
#' @param microbeNames vector of strings of microbial group names
#' @return group name (string) if xname is a strain name. If xname is not a the name of a strain it will simply return xname unchanged. 
#' @export



getGroupName = function(xname, microbeNames) {
    
    g = as.numeric(gregexpr(pattern = "\\.", xname)[[1]]) - 1  #find where dot(s) is(are)
    # gsub('^([^\\.]+)\\.?.*$','\\1', xname) (alternative method)
    if (max(g) <= 0) {
        # not a strain name, therefore do not alter xname
        x = xname
    } else {
        gname = substring(xname, first=1, last=max(g))
        if (gname %in% microbeNames) {
            x = gname
        } else {
            x = xname
        }
    }
    
    
    return(x)
    
}
