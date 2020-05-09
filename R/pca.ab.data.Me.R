#' A superfunction to estimate correlations and relationships among variables in a data set. 
#'
#' This function estimates a Spearman's correlation matrix, an abs(distance matrix), a hclust tree, a PCA from the correlation matrix, and numerous estimates of the effective number of variables in a data set (Me)
#' @param variabledata a vector of class categorical, for all elements tested, identifying a category | or bin | box that each element belong to
#' @keywords correaltion matrix, distance matrix, hclust, Me
#' @export
#' @examples
#' pca.ab.data.Me()
pca.ab.data.Me = function( wdata){
    cat("- Estimating principal components from the abundance data\n")
    ## number of variables
    M = ncol(wdata)

    ## estimate missingness
    cat("\t1. Estimating feature missingness\n")
    ##
    fmis = apply(wdata, 2, function(x){ sum(is.na(x))/length(x) })

    ## filter out features with missingness >5%
    cat("\t2. Removing feature with missingness > 5%\n")
    ##
    w = which(fmis > 0.05)
    if(length(w)> 0){
        wdata = wdata[, -w]
        Mf = ncol(wdata)
    } else {
        Mf = M
    }

    ## impute missingness to the median
    cat("\t3. Imputing missingness to the median\n")
    ##
    for(i in 1:ncol(wdata)){
        x = wdata[, i]
        m = median(x, na.rm = TRUE)
        x[is.na(x)] = m
        wdata[, i] = x
    }

    ## z-transform the data
    cat("\t4. z-transform the data\n")
    wdata = apply(wdata, 2, ztransform)

    ## PCA
    cat("\t5. estimating principal components from correlation matrix\n")
    pca = prcomp(wdata, center = FALSE, scale = FALSE)
        #eigenvalues = summary(pca)[[6]][2, ]
    
    cat("\t6. extract eigenvalues / variance explained by PCs\n")
    varexp = t(summary(pca)[[6]][2:3, ])
    eigenvalues = varexp[,1]

    ## effective number of markers
    cat("\t7. Estiamte Me (the effective number of markers)\n")
        ## Goa 2008
    m95 = length(which(varexp[,2] < 0.95))+1
    m995 = length(which(varexp[,2] < 0.995))+1
    
        ## Li and Ji 2005
    ve = eigenvalues
    w = which(ve > 0.01)
    ve[w] = 1
    LJ_Me = round( sum(ve) )
    ##
    simpleM = data.frame( M = M, LiJi_Me = LJ_Me,  Goa_M_95 = m95  , Goa_M_995 = m995  )
    
    ## data out
    out = list(simpleM = simpleM, 
        varexp = varexp,
        pca = pca,
        variabledata = wdata )
    ## return data
    return(out)
}