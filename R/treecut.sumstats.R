#' A function to estimate some intra-cluster summary statstics in you hclust object.
#'
#' This function estimates a few summary statistics (N, min, mean, max) for each cluster, with an N > 1, at your defined cutheight in your hclust object.
#' @param tree a hclust object of you variables, built with the cormat and the hclust method "complete"
#' @param cormat a correlation matrix among you variables
#' @param cutheight the height at which you wish to cut your tree 
#' @keywords tree cut, summary statistics
#' @export
#' @examples
#' treecut.sumstats()
treecut.sumstats = function( tree, cormat, cutheight ){
  ##############  check for hclust object
  if(class(tree) != "hclust"){
    stop("You must pass a hclust object as the tree variable")
  }
  ##############  check for a correlation matrix
  if( class(cormat) != "matrix"  | length(tree$labels) != nrow(cormat) ){
    stop( paste0("You must pass the correlation matrix - NOT distance matrix - used to generate your hclust object.\nThe number of terminal branches must equal the number of rows and columns of matrix") )
  }
  ##############  check for a single numeric float
  if( length(cutheight) != 1  ){
    stop(  "You must pass a single cutheight float >= 0  & <= 1." )
  }
  if( class(cutheight) != "numeric" |  cutheight < 0 | cutheight > 1 ){
    stop("You must pass a single, cutheigth float no smaller than 0 and no larger than 1.")
  }
  ############## Produce a warning if the hclust method used was not "complete".
  if( !tree$method %in% c( "complete" ,"average","mcquitty") ){
    warning( paste0("It is advised to use the hclust method 'complete', 'average', or 'mcquitty.' You used the method ", tree$method, ".") )
  }
  ############## perform tree cut
  k = stats::cutree(tree, h = cutheight)
  ############## unique K clusters with more than one variable in them
  uniqueKs = names( which( table(k) > 1 ) )
  ##############  iterate over these clusters and estimate sumstats
  sumstats = t(
    sapply(uniqueKs, function(i){
        i = as.numeric(i)
        w = which(k == i)
        n = names(k[w])
        tempd = as.dist(abs(cormat[n,n]))
        ss_out = data.frame( k = i, N = length(w), min_rho = min(tempd, na.rm = TRUE) , 
          mean_rho =  mean(tempd, na.rm = TRUE) , max_rho = max(tempd, na.rm = TRUE), sd_rho = sd(tempd, na.rm = TRUE) )
        return(unlist(ss_out))
    } ) ## end of sapply
  )  ## end of t() function
  out = list( sumstats = sumstats, k = k) 
  return(out)

}
