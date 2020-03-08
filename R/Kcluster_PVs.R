#' A function to identify the Principal Variable of a group
#'
#' A function to identify the Principal Variable of a group. The final step in the iPVs pipeling to identify the final PV to stand as the tagging and independent variable for the cluster of variables. 
#' @param variabledata a data.frame of the variable data used to build your cormat and hclust tree
#' @param Kmembers a list, where each object in the list is a vector of variable IDs, identifying the members of a cluster or group.
#' @param myPVs a vector of my Principal Variable (myPVs) for each cluster
#' @keywords principal variables, tree cut, iPV
#' @export
#' @examples
#' Kcluster_PVs()
Kcluster_PVs = function( variabledata, Kmembers, myPVs ){
  ## summary table of top PV for each cluster
  PVtable = c()
  ## list for all dat
  PVresults = list()
  
  for(i in 1:length(Kmembers)){
    print(i)
    ## variable IDs
    n = Kmembers[[i]]
    ###
    if(length(n) > 1){
      ## temporary data
      tempd = variabledata[, n]
      ## PVA analysis 
      PV = VarRep(tempd)
      w = which( PV$variable == myPVs[i] )
      PV2pass = data.frame( PV[w, c(1,3,4) ], PVArank = w )
      ###
      PVtable = rbind( PVtable, PV2pass )
    } else {
      PV = data.frame(variable = n, initial_sumR2 = 1, VarExp_individually = 1, added_vexp = 1,  cum_vexp = 1 )
      PV2pass = data.frame( PV[, c(1,3,4) ], PVArank = 1 )
      ##
      PVtable = rbind( PVtable, PV2pass )
      ##
    }
    ####
    PVresults[[i]] = PV
  }
  ## add names to rows and lists
  rownames(PVtable) = myPVs
  names(PVresults) = myPVs
  ##
  out = list(PVtable = PVtable, PVresults = PVresults)
  return(out)
}


