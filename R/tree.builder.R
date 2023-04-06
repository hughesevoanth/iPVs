#' A superfunction to estimate correlations and relationships among variables in a data set. 
#'
#' This function estimates a Spearman's correlation matrix, an abs(distance matrix), a hclust tree, a PCA from the correlation matrix, and numerous estimates of the effective number of variables in a data set (Me)
#' @param variabledata a vector of class categorical, for all elements tested, identifying a category | or bin | box that each element belong to
#' @keywords correaltion matrix, distance matrix, hclust, Me
#' @export
#' @examples
#' tree.builder()
tree.builder = function( wdata, cor_method = "spearman", dist_method = "R", hclust_meth = "average", cmat = NULL ){
    
    ## number of variables
    M = ncol(wdata)

    ## estimate correlation matrix
    if( is.null(cmat) ){
      ## number of variables
      M = ncol(wdata)
      
      cat("1. estimating pairwsie correlations\n")
        cmat = cor(wdata, use = "pairwise.complete.obs", method = cor_method )
        ## turn NAs into 0 values
        cmat[is.na(cmat)] = 0
    } else {
      cat("1. pairwsie correlation matrix provided by user\n")
      M = ncol(cmat)
    }

    ## estimate a distance matrix
    cat("2. constructing distance matrix\n")
    if(dist_method == "R2"){
        dmat = as.dist(1 - (cmat * cmat) )    
    } else {
        if(dist_method == "R"){
            dmat = as.dist(1 - abs( cmat ) )    
        } else { 
            stop("The dist_method must be either 'R' (1-abs(R)) or 'R2' (R*R).")
        }
    }
    

    ## generate a hclust dendrogram
    cat("3. generating dendrogram\n")
    if(hclust_meth == "complete"){
        tree = hclust(dmat,  method = "complete")
    } else {
        if(hclust_meth == "average"){
            tree = hclust(dmat,  method = "average")
        } else { 
            if(hclust_meth == "mcquitty"){
                tree = hclust(dmat,  method = "mcquitty")
            } else { 
                stop("The hclust_meth must be either 'complete' (max dist), 'average' (UPGMA), or 'mcquitty' (WPGMA).")
            }
        }
    }

    ## eigenvalues
    cat("4. estimating eigenvalues from correlation matrix\n")
    eigenvalues = eigen(cmat, symmetric = TRUE)$values

    ## PCA
    cat("5. estimating principal components from correlation matrix\n")
    pca = prcomp(cmat, center = TRUE, scale = TRUE)
        #eigenvalues = summary(pca)[[6]][2, ]
    varexp = t(summary(pca)[[6]][2:3, ])
    
    ## effective number of markers
        ## Goa 2008
    m95 = length(which(varexp[,2] < 0.95))+1
    m995 = length(which(varexp[,2] < 0.995))+1
        ## Cheverud 2001
    eigen_var = var(eigenvalues)
    Meff = round( 1 + (M-1) * (1 - (eigen_var / M)) )
        ## Li and Ji 2005
    ve = varexp[,1]
    w = which(ve > 0.01)
    ve[w] = 1
    LJ_Me = round( sum(ve) )
    ##
    simpleM = data.frame( Cheverud_Me = Meff, LiJi_Me = LJ_Me,   Goa_M_95 = m95 , Goa_M_995 = m995   )
    
    ## data out
    out = list( variabledata = wdata, 
        cormat = cmat, distmat = dmat, tree = tree, 
        pca = pca, eigenvalues = eigenvalues,  varexp = varexp, simpleM = simpleM
       )
    ## return data
    return(out)
}


