#' A function to identify and extract the principal variables in your dataset. 
#'
#' A function to identify and extract the principal variables in your dataset. The function will iterate over each newly generated list of PVs to insure that they do meet your tree cutheight standards.
#' @param variabledata a data.frame of the variable data used to build your cormat and hclust tree
#' @param tree a hclust object of you variables, built with the cormat and the hclust method "complete"
#' @param cormat your correlation matrix
#' @param distmat your distance matrix
#' @param cutheight the height at which you wish to cut your tree 
#' @keywords principal variables, tree cut
#' @export
#' @examples
#' ind.pvs()
ind.pvs = function( variabledata, tree, cormat, distmat, cutheight){
  TC_iteration = list()
  nomoreclusters = 0
  ### while loop
  ### while there are still clusters, at height h, in each newly generated tree using the PVs
  ### thencut tree and make new PVs again. 
  i = 1
  while(nomoreclusters == 0){
    
    cat("Tree cut and PV identifier iteration ", i, "\n")
    ## set a seed
    set.seed(20200224)
    ## identify the Principal Variables (PVs)
    PVs = treecut.pvs( tree = tree, variabledata = variabledata, cutheight = cutheight)
    
    ## add the k cluster identifiers to the list object Kclusters
    title = paste0("treecut.pvs_iteration_", i)
    TC_iteration[[title]] = PVs

    ## vector of Principal variable names
    n = as.character(PVs$pvs[, 3])
    
    ## create a new temporary distance matrix using the Newly identified PVs
    newdistmat = as.dist(as.matrix( distmat )[n, n])
    
    ## estiamte a new hclust tree
    tree = hclust(newdistmat,  method = "complete")
    
    ## estimate the number of clusters in this new tree using only the PVs
    k = table( stats::cutree(tree, h = cutheight) )
    ###############
    ## if the number of clusters in the tree cut equals the number of variables then we are done. 
    ## kill the while loop
    if( length(n) == length( k ) ){
      nomoreclusters = 1
    } else {
      i = i+1  
    }
    ###############
    
  
  }
  #### estimate the maximum correlation that remains in the data set
  n = as.character(PVs$pvs[, 3])
  MaxCor_remaining = max( as.dist( abs( cormat[n,n] ) ) , na.rm = TRUE)
  ###
  out = list(PVs = PVs, tree = tree, cutheight = cutheight, MaxCor_remaining = MaxCor_remaining, treecut_iterations = TC_iteration)
  return(out)
}