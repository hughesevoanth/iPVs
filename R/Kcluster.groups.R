#' A function to identify all members of a cluster, generated from iterative tree cuts
#'
#' A function to identify all members of a cluster, generated from iterative tree cuts
#' @param ks a list of k ( cutree() ) group identifiers for each iterative tree cut performed to generate a fully independent listof PVs.
#' @param pvs a vector of identified PVs (principal variables)
#' @keywords principal variables, tree cut, independent variables
#' @export
#' @examples
#' Kcluster.groups()
Kcluster.groups = function( ind_pv_iterations ){
  final_iPV = as.character( ind_pv_iterations[[length(ind_pv_iterations)]]$pvs$PV )

  ###########################
  ##  find cluster members 
  ## for each Principal Variable
  ############################
  cluster_members = lapply(final_iPV, function(pv){
    # print(pv)
    
    ## a new vector to store
    ## cluster members, starting with 
    ## the (current) iPV
    k_members = pv
    
    ### LOOP OVER tree cut iterations
    ### to find iterative cluster members
    for(i in length(ind_pv_iterations):1 ){
      ## k data for iteration "i"
      k = ind_pv_iterations[[i]]$k
      
      ## this|these cluster-members belong to which cluster(s) ?
      w = which( names( k ) %in% k_members )
      k_ids = k[w]
      ## identify all variables with this k_id
      w = which( k %in% k_ids )
      new_k_members = as.character( names(k)[w] )
      ## add those members to the list (this includes itself)
      k_members = unique( c(k_members, new_k_members) )  
    }
    ### END LOOP

    ## return
    return( k_members )

    }) ## end of lapply
}
