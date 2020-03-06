#' A super function that runs through all of the necessary staeps for you to identify Independent Principal Variables
#'
#' A super function that runs through all of the necessary staeps for you to identify Independent Principal Variables
#' @param variabledata a data.frame of the variable data used to build your cormat and hclust tree
#' @param cutheight tree cut height. Value is equal to a dissimilarity of 1 - Spearman's rho.
#' @keywords data reduction, independent variables, principal variables, tree cut
#' @export
#' @examples
#' iPVs()
iPVs = function( variabledata, cutheight = 0.5 ){
  ## estiamte correlation matrix, build tree, generate PCA from correlation matrix
  cat(paste0("(I) tree.builder -- \n"))
  wdata = tree.builder(variabledata)
  
  ## identify the PVs (independent principal variables)
  cat(paste0("(II) ind.pvs -- identify independent clusters and initial principal variables\n"))
  StudyPVs = ind.pvs( variabledata = wdata$variabledata,
  tree = wdata$tree,
  cormat = wdata$cormat,
  distmat = wdata$distmat,
  cutheight = cutheight )
  
  initial_ind_PVs = as.character( StudyPVs$PVs$pvs[,"PV"] )

  ## identify all variables that belong to a group. i.e. each PV is tagging which other variables? 
  cat(paste0("(III) Kcluster.groups -- identify all variables|members of a cluster, iteratively.\n"))
  PV_cluster_members = Kcluster.groups( ind_pv_iterations = StudyPVs$treecut_iterations )
  
  ## re-estiamte PV and the VarExp by that top PV for the total variation of group members.
  cat(paste0("(IV) Kcluster_PVs -- identify the final set of PV for each cluster, and estimate the variance explained.\n"))
  NewPV = Kcluster_PVs(variabledata = wdata$variabledata, Kmembers = PV_cluster_members )
  final_ind_PVs = as.character(NewPV$PVtable[,1])
  vexp = NewPV$PVtable[,2]
  
  ## place all of the useful data into a table
  cat(paste0("(V) Generate summary table with PVs, cluster members, and variance explained.\n"))
  clustersize = unlist( lapply(PV_cluster_members, length) )
  groupmembers = unlist( lapply(PV_cluster_members, function(x){ paste(x, collapse = ":") } ) )
  
  iPV_table = data.frame(PVs = final_ind_PVs , 
    clustersize = clustersize,
    VarExp_by_PV = vexp,
    groupmembers = groupmembers )

  ### data out
  out = list(iPV_table = iPV_table, 
    PV_cluster_members = PV_cluster_members,
    PVresults = NewPV$PVresults,
    workingdata = wdata )
  return(out)

}