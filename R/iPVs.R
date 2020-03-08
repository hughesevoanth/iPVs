#' A super function that runs through all of the necessary staeps for you to identify Independent Principal Variables
#'
#' A super function that runs through all of the necessary staeps for you to identify Independent Principal Variables
#' @param variabledata a data.frame of the variable data used to build your cormat and hclust tree
#' @param cor_method the correlation method used to construct the correlation matrix
#' @param dist_method the distance method to generate the distance matrix
#' @param hclust_meth the clustering method used by hclust
#' @param cutheight tree cut height. Value is equal to a dissimilarity of 1 - Spearman's rho.
#' @keywords data reduction, independent variables, principal variables, tree cut
#' @export
#' @examples
#' iPVs()
iPVs = function( variabledata, cor_method = "spearman", dist_method = "R2", hclust_meth = "average", cutheight = 0.5 ){

  ## estiamte correlation matrix, build tree, generate PCA from correlation matrix
  cat(paste0("(I) tree.builder -- \n"))
  wdata = tree.builder(variabledata, cor_method = cor_method, dist_method = dist_method, hclust_meth = hclust_meth )
  
  PVlistout = lapply(cutheight, function(CH){    
      ## identify the PVs (independent principal variables)
      cat(paste0("(II) ind.pvs -- identify independent clusters and initial principal variables\n"))
      StudyPVs = ind.pvs( variabledata = wdata$variabledata,
      tree = wdata$tree,
      cormat = wdata$cormat,
      distmat = wdata$distmat,
      cutheight = CH, 
      hclust_meth = hclust_meth )

      ## the last tree cut iteration defines our 
      ## data sets PVs
      ## the PVs identified for your data set
      mypvs = as.character( StudyPVs$PVs$pvs[,"PV"] )

      ## identify all variables that were clustered with your PV in any one of the tree cut iteractions
      ##  i.e. your PVs are tagging which other variables in your data set?
      cat(paste0("(III) Kcluster.groups -- identify all variables|members of a cluster, iteratively.\n"))
      PV_cluster_members = Kcluster.groups( ind_pv_iterations = StudyPVs$treecut_iterations )
      
      ## Perform a  PVA for the last time for for each of your iterative-super-clusters
      ##  and the VarExp by that top PV for the total variation of group members.
      cat(paste0("(IV) Kcluster_PVs -- identify the final set of PV for each cluster, and estimate the variance explained.\n"))
      Final_PVA_results = Kcluster_PVs( variabledata = wdata$variabledata, Kmembers = PV_cluster_members, myPVs = mypvs )
      
      ##############################
      ## place all of the useful data
      ## in a final table
      ##############################
      #cat(paste0("(V) Generate summary table with PVs, cluster members, and variance explained.\n"))
      #clustersize = unlist( lapply(PV_cluster_members, length) )
      #groupmembers = unlist( lapply(PV_cluster_members, function(x){ paste(x, collapse = ":") } ) )
      
      cat(paste0("(V) Generate summary table with PVs, and variance explained.\n"))
      clustersize = unlist( lapply(Final_PVA_results$PVresults, nrow) )
      
      
      iPV_table = data.frame(PVs = Final_PVA_results$PVtable$variable , 
        clustersize = clustersize,
        VarExp_by_PV = Final_PVA_results$PVtable$VarExp_individually,
        PVArank = Final_PVA_results$PVtable$PVArank )
        #groupmembers = groupmembers )

      ### data out
      out = list(iPV_table = iPV_table, 
        PV_cluster_members = PV_cluster_members,
        PVresults = Final_PVA_results$PVresults,
        workingdata = wdata )

      return(out)
    })
  ## add names to each list
  names(PVlistout) = paste0("cutheight_", cutheight)

  ### reduce if there is only 1 cut height
  if(length(PVlistout) == 1){
    PVlistout = PVlistout[[1]]
  }

  ##
  return(PVlistout)
}