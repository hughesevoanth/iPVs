#' A function to iterate over various tree cut height and extract the minimum observed correaltion coefficient within a cluster.
#'
#'  A function to iterate over various tree cut height and extract the minimum observed correaltion coefficient within a cluster.
#' @param tree a hclust object of you variables, built with the cormat and the hclust method "complete"
#' @param cormat your data.frame of variables used to generate your tree
#' @param cutheights the heights at which you wish to cut your tree 
#' @keywords treecut.sumstats, minimum rho or correlation coefficient
#' @export
#' @examples
#' treecut.iterator.4minrho()
treecut.iterator.4minrho = function( tree, cormat, cutheights){
  ## sapply function
  cutheight_MinCor = t( sapply(cutheights, function(h){
    
    ## summary statistics for each specific tree cut height
    ss = treecut.sumstats(tree = tree, cormat = cormat, cutheight = h)
    
    ## define the minimum observed intra-group rho at this tree cut height
    if( length( table( ss$k ) ) != length(ss$k) ){
      obs_min = signif( min( ss$sumstats[, 3], na.rm = TRUE) , d = 4)
      } else { 
        obs_min = 1 
      }
    
    ## values to retrun
    out = c(h, obs_min)
    names(out) = c("cutheight","obs_minimum_rho")
    return(out)  
    
    }) )

  ###
  return(cutheight_MinCor)
}

