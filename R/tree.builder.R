#' A superfunction to estimate correlations and relationships among variables in a data set. 
#'
#' This function estimates a Spearman's correlation matrix, an abs(distance matrix), a hclust tree, a PCA from the correlation matrix, and numerous estimates of the effective number of variables in a data set (Me)
#' @param variabledata a vector of class categorical, for all elements tested, identifying a category | or bin | box that each element belong to
#' @keywords correaltion matrix, distance matrix, hclust, Me
#' @export
#' @examples
#' tree.builder()
tree.builder = function( wdata ){
    
    ## number of variables
    M = ncol(wdata)

    ## estimate correlation matrix
    cat("1. estimating pairwsie correlations\n")
    cmat = cor(wdata, use = "pairwise.complete.obs", method = "spearman")

    ## estiamte a distance matrix
    cat("2. constructing distance matrix\n")
    dmat = as.dist(1 - abs( cmat ) )

    ## generate a hclust dendrogram
    cat("3. generating dendrogram\n")
    tree = hclust(dmat,  method = "complete")

    ## eigenvalues
    cat("4. estimating eigenvalues\n")
    eigenvalues = eigen(cmat, symmetric = TRUE)$values

    ## PCA
    cat("5. estimating principal components\n")
    pca = prcomp(cmat, center = TRUE, scale = TRUE)
    
    varexp = t(summary(pca)[[6]][2:3, ])
    
    m95 = length(which(varexp[,2] < 0.95))+1
    m995 = length(which(varexp[,2] < 0.995))+1
    Me = M*(1 - ( M-1 ) * ( var(eigenvalues)/M^2 ) )
    ##
    simpleM = data.frame( M_95 = m95, M_995 = m995, Me = Me )
    
    ## data out
    out = list(variabledata = wdata, cormat = cmat, distmat = dmat, tree = tree, 
      eigenvalues = eigenvalues, pca = pca, varexp = varexp, simpleM = simpleM)
    ## return data
    return(out)
}
