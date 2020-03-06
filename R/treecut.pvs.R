#' A function to identify the principal variables in your variable data set, given a defined tree cut height. 
#'
#' This function to identify the principal variables in your variable data set, given a defined tree cut height. 
#' @param tree a hclust object of you variables, built with the cormat and the hclust method "complete"
#' @param variabledata your data.frame of variables used to generate your tree
#' @param cutheight the height at which you wish to cut your tree 
#' @keywords tree cut, principal variables
#' @export
#' @examples
#' treecut.pvs()
treecut.pvs = function( tree, variabledata , cutheight){
  ##############  check for hclust object
  if(class(tree) != "hclust"){
    stop("You must pass a hclust object as the tree variable")
  }
  ##############  check for presence of all tree labels in variabledata
  if( sum(!tree$labels %in%  colnames(variabledata) ) ){
    stop( paste0("Make sure that all of the variables in your hclust object are also in your variabledata data.frame.") )
  }
  ##############  check for a single numeric float
  if( length(cutheight) != 1  ){
    stop(  "You must pass a single cutheight float >= 0  & <= 1." )
  }
  if( class(cutheight) != "numeric" |  cutheight < 0 | cutheight > 1 ){
    stop("You must pass a single, cutheigth float no smaller than 0 and no larger than 1.")
  }
  ############## Produce a warning if the hclust method used was not "complete".
  if( tree$method != "complete" ){
    warning( paste0("It is advised to use the hclust method complete. You used the method ", tree$method, ".") )
  }
  ###############
  ## cut tree
  k = stats::cutree(tree, h = cutheight)
  ## vector of unique Ks
  singleKs = as.numeric( names( which( table(k) == 1 ) ) )
  ## vector of unique Ks with n > 1
  groupKs = as.numeric( names( which( table(k) > 1 ) ) )
  ## iterate of goupKs to identify a Principal Variable
  if(length(groupKs) > 0){
  pvs = t( sapply(groupKs, function(i){
    n = names( which(k == i) )
    tempd = variabledata[, n]
    ##
    topPV =  PVA( names(tempd), tempd)
    out = data.frame(k = i, clusN = table(k)[i],  PV = as.character(topPV$Selected)[1], VarExp = topPV$Added.Propn[1])
    return(out)
    }) )
  ## redefine pvs as a proper data.frame without built in lists :-( 
  pvs = data.frame(k = as.numeric(unlist(pvs[,1])), clusN = as.numeric(unlist(pvs[,2])), 
    PV = as.character(unlist(pvs[,3])), VarExp =  as.numeric(unlist(pvs[,4])))
  } else {
    pvs = data.frame( )
  }
  ##  build data frame for all single variable K groups
  w = which(k %in% singleKs)
  n = names(k)[w]
  singleK_pvs = data.frame(k = singleKs, clusN = 1, PV = n, VarExp = 1)
  ## combine the pvs
  pvs = rbind(pvs, singleK_pvs)
  ## order the pvs by k
  o = order(pvs[,1])
  pvs = pvs[o, ]
  rownames(pvs) = paste0("K", pvs[,1])
  ##
  return( list( pvs = pvs, k = k ) )
}
