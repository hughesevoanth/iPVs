#' A function to identify the Principal Variable of a group
#'
#' A function to identify the Principal Variable of a group. The final step in the iPVs pipeling to identify the final PV to stand as the tagging and independent variable for the cluster of variables. 
#' @param variabledata a data.frame of the variable data used to build your cormat and hclust tree
#' @param Kmembers a list, where each object in the list is a vector of variable IDs, identifying the members of a cluster or group.
#' @keywords principal variables, tree cut, iPV
#' @export
#' @examples
#' Kcluster_PVs()
Kcluster_PVs = function( variabledata, Kmembers ){
  ## summary table of top PV for each cluster
  PVtable = c()
  ## list for all dat
  PVresults = list()
  
  for(i in 1:length(Kmembers)){
    ## variable IDs
    n = Kmembers[[i]]
    if(length(n) > 1){
      ## temporary data
      tempd = variabledata[, n]
      ## PVA analysis
      PV =  PVA( names(tempd), tempd)
      PV$Selected = as.character(PV$Selected)  
    } else {
      PV = data.frame(Variable = 1, Selected = n, h.partial = NA, Added.Propn = 1, Cumulative.Propn = 1)
      #names(PV) = c("Variable","Selected","h.partial","Added.Propn","Cumulative.Propn")
    }
    ####
    PVtable = rbind( PVtable, PV[1, c(2,4)] )
    PVresults[[i]] = PV
  }
  out = list(PVtable = PVtable, PVresults = PVresults)
  return(out)

}