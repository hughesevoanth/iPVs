#' A function to identify the Principal Variable of a group
#'
#' A function to identify the Principal Variable of a group. The final step in the iPVs pipeling to identify the final PV to stand as the tagging and independent variable for the cluster of variables. 

#' @param cmat a correlation matrix you can pass to the function, which will skip the generation of a new matrix.
#' @param Kmembers a list, where each object in the list is a vector of variable IDs, identifying the members of a cluster or group.
#' @param myPVs a vector of my Principal Variable (myPVs) for each cluster
#' @keywords principal variables, tree cut, iPV
#' @export
#' @examples
#' Kcluster_PVs()
Kcluster_PVs = function( cmat, Kmembers, myPVs ){
  ## summary table of top PV for each cluster
  PVtable = c()
  ## list for all dat
  PVresults = list()
  
  for(i in 1:length(Kmembers)){
    # print(i)
    ## variable IDs
    n = Kmembers[[i]]
    ###
    if(length(n) > 1){
      ## temporary data
      tempd = cmat[n, n]
      ## PVA analysis 
      PV = VarRep(tempd)
      ## Which PV is my PV ?
      w = which( PV$variable == myPVs[i] )
      ## rank by Variance Explained individually
      o = order( PV$VarExp_individually, decreasing = TRUE )
      PV = PV[o,]
      ## Redefine what the PV is.
      myPVs[i] = PV$variable[1]
      # temp = PV[o,]
      # q = which(temp$variable == myPVs[i] ); rm(temp)
      ##
      # PV2pass = data.frame( PV[w, c(1,3,4) ], PVArank = w, ind_VarExp_rank = q  )
      PV2pass = data.frame( PV[1, c(1,3,4) ], PVArank = 1, ind_VarExp_rank = 1  )
      ###
      PVtable = rbind( PVtable, PV2pass )
    } else {
      PV = data.frame(variable = n, initial_sumR2 = 1, VarExp_individually = 1, added_vexp = 1,  cum_vexp = 1 )
      PV2pass = data.frame( PV[, c(1,3,4) ], PVArank = 1, ind_VarExp_rank = 1 )
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


